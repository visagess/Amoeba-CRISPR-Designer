"""Small local form; all analysis runs in a separate process."""
from runtime import storage_path, settings
import os
from pathlib import Path
import queue
import subprocess
import sys
import threading
import json
import datetime as dt
import re
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

ROOT = Path(__file__).resolve().parent


def main():
    window = tk.Tk()
    window.title("Amoeba CRISPR Designer · Cas12a RNA 设计")
    window.geometry("1000x980")
    frame = ttk.Frame(window, padding=24)
    frame.pack(fill="both", expand=True)
    ttk.Label(frame, text="输入名称，使用 CaSilico 和 EasyDesign 设计 Cas12a crRNA", font=("Microsoft YaHei UI", 16)).pack(anchor="w")
    ttk.Label(frame, text="目标种 + 背景属 → 候选特异/保守区域 → CaSilico → RNA 和结果报告").pack(anchor="w", pady=(4, 18))
    form = ttk.Frame(frame)
    form.pack(fill="x")
    name, database, groups, threads = tk.StringVar(value="Naegleria gruberi"), tk.StringVar(value="NCBI"), tk.StringVar(value="10"), tk.StringVar(value="8")
    mode=tk.StringVar(value="18S在线搜索")
    genera=tk.StringVar(value="Entamoeba, Naegleria")
    background = tk.StringVar()
    mismatches=tk.StringVar(value="2")
    coverage_min=tk.StringVar(value='')
    coverage_max=tk.StringVar(value='')
    local_inputs=[]
    background_min=tk.StringVar(value='2')
    max_exceptions=tk.StringVar(value='3')
    min_bp=tk.StringVar(value='800')
    max_bp=tk.StringVar(value='4000')
    fields=[("分析方式",mode),("目标种拉丁学名",name),("数据库",database),("背景属（英文逗号分隔）",genera),
            ("18S最短长度（bp，含边界）",min_bp),("18S最长长度（bp，含边界）",max_bp),
            ("18S窗口允许错配数（0–4）",mismatches),("保守覆盖物种下限（空=全部）",coverage_min),("保守覆盖物种上限（空=全部）",coverage_max),
            ("特异区与背景最少错配数（含边界）",background_min),("每段允许不特异的背景物种数",max_exceptions),
            ("CaSilico 设计区域数",groups),("线程数",threads),("额外背景 FASTA（可选）",background)]
    for row, (label, variable) in enumerate(fields):
        ttk.Label(form, text=label).grid(row=row, column=0, sticky="w", padx=(0, 16), pady=5)
        if variable is database:
            entry=ttk.Combobox(form,textvariable=variable,state="readonly",values=["NCBI","AmoebaDB（公共存档）"])
        elif variable is mode:
            entry=ttk.Combobox(form,textvariable=variable,state="readonly",values=["18S在线搜索","18S本地FASTA","特异区域 + CaSilico","共有基因 + CaSilico"])
        else:
            entry=ttk.Entry(form,textvariable=variable)
        entry.grid(row=row, column=1, sticky="ew", pady=5)
    form.columnconfigure(1, weight=1)
    ttk.Button(form, text="选择", command=lambda: background.set(filedialog.askopenfilename(filetypes=[("FASTA", "*.fa *.fasta *.fna *.gz"), ("全部", "*")]))).grid(row=len(fields)-1, column=2, padx=(8, 0))
    ttk.Label(frame,text="在线模式可自定目标种、背景属和18S长度。目标属自动纳入，用于属内保守比较及背景筛查。\n长度范围按提取后的核18S计算，适用于全部目标、背景及后续设计。本地输入每个文件一个物种。\n覆盖上下限留空=本轮目标属全部已定种；四组结果含GC及物种名单，EasyDesign独立预测21nt候选。",wraplength=880).pack(anchor="w",pady=(12,10))
    buttons = ttk.Frame(frame)
    buttons.pack(fill="x", pady=8)
    text = tk.Text(frame, wrap="word", state="disabled", font=("Consolas", 10))
    text.pack(fill="both", expand=True, pady=(8, 0))
    events = queue.Queue()
    state = {"process": None, "output": None}

    def append(line):
        text.configure(state="normal")
        text.insert("end", line)
        text.see("end")
        text.configure(state="disabled")

    def start():
        try:
            n, t = int(groups.get()), int(threads.get())
            if n < 0 or t < 1 or not name.get().strip():
                raise ValueError()
        except ValueError:
            messagebox.showerror("请检查输入", "请填写物种学名；区域数至少为 0，线程数至少为 1。")
            return
        command=[sys.executable.replace("pythonw.exe","python.exe"),"-u"]
        db="ncbi" if database.get()=="NCBI" else "amoebadb"
        if mode.get() in ("18S在线搜索","18S本地FASTA"):
            supplied=[s.strip() for s in genera.get().replace('，',',').split(',') if s.strip()]
            if background.get():
                messagebox.showerror("18S本地输入","18S文件请通过“导入FASTA并选目标文件”确认每个物种及目标文件。")
                return
            try:
                from ssu18s_policy import length_limits
                minimum_bp,maximum_bp=length_limits(int(min_bp.get()),int(max_bp.get()))
                mismatch_limit=int(mismatches.get())
                coverage_lower=int(coverage_min.get()) if coverage_min.get().strip() else None
                coverage_upper=int(coverage_max.get()) if coverage_max.get().strip() else None
                background_distance=int(background_min.get())
                exception_limit=int(max_exceptions.get())
                if not 0<=mismatch_limit<=4:raise ValueError()
                if any(v is not None and v<1 for v in (coverage_lower,coverage_upper)):raise ValueError()
                if coverage_lower is not None and coverage_upper is not None and coverage_lower>coverage_upper:raise ValueError()
                if not 1<=background_distance<=20:raise ValueError()
                if exception_limit<0:raise ValueError()
            except ValueError:
                messagebox.showerror('请检查筛选参数','20 ≤ 最短bp ≤ 最长bp；目标错配0–4；覆盖下限≥1且不大于上限；背景最少错配1–20；例外种数≥0。')
                return
            command.extend([str(ROOT/'ssu18s_run.py'),'--threads',str(t),'--mismatches',str(mismatch_limit),
                            '--background-min-mismatches',str(background_distance),'--max-exception-species',str(exception_limit),
                            '--min-bp',str(minimum_bp),'--max-bp',str(maximum_bp)])
            if mode.get()=="18S在线搜索":
                command.extend(['--target',name.get().strip(),'--database',db])
                for genus in supplied:command.extend(['--background-genus',genus])
            if coverage_lower is not None:command.extend(['--coverage-min',str(coverage_lower)])
            if coverage_upper is not None:command.extend(['--coverage-max',str(coverage_upper)])
            if mode.get()=="18S本地FASTA":
                if not local_inputs:
                    messagebox.showerror('需要本地FASTA','请先点击“导入FASTA并选目标文件”。');return
                folder=Path(str(storage_path('local_inputs')));folder.mkdir(parents=True,exist_ok=True)
                manifest=folder/(dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f')+'.json')
                manifest.write_text(json.dumps({'files':local_inputs},ensure_ascii=False),encoding='utf-8')
                command.extend(['--local-manifest',str(manifest)])
        elif mode.get()=="特异区域 + CaSilico":
            genus_names=[s.strip() for s in genera.get().replace('，',',').split(',') if s.strip()]
            if not genus_names:
                messagebox.showerror("需要背景范围","请填写至少一个背景属的拉丁名称。")
                return
            command.extend([str(ROOT/'region_pipeline.py'),'--target',name.get().strip(),'--database',db,
                            '--design-regions',str(n),'--threads',str(t)])
            for genus in genus_names:
                command.extend(['--background-genus',genus])
            if background.get():command.extend(['--background-fasta',background.get()])
        else:
            command.extend([str(ROOT/'amoebacrisp.py'),'--name',name.get().strip(),'--database',db,
                            '--max-groups',str(n),'--threads',str(t)])
            if background.get():command.extend(['--background',background.get()])
        launch.configure(state="disabled")
        append("\n开始分析……\n")

        def work():
            try:
                env = dict(os.environ, PYTHONUTF8="1", PYTHONIOENCODING="utf-8")
                process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace",
                                           creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
                state["process"] = process
                for line in process.stdout:
                    events.put(("line", line))
                events.put(("done", process.wait()))
            except Exception as exc:
                events.put(("line", str(exc) + "\n"))
                events.put(("done", 1))
        threading.Thread(target=work, daemon=True).start()

    def poll():
        try:
            while True:
                kind, value = events.get_nowait()
                if kind == "line":
                    append(value)
                    if value.startswith("输出目录："):
                        state["output"] = value.split("：", 1)[1].strip()
                else:
                    launch.configure(state="normal")
                    state["process"] = None
                    append("\n分析完成，请打开结果文件夹。\n" if value == 0 else "\n分析未完成，请查看上方错误及输出文件夹中的日志。\n")
        except queue.Empty:
            pass
        window.after(150, poll)

    def open_output():
        folder = state["output"] or (str(storage_path('outputs')) if mode.get().startswith('18S') else
                                    str(storage_path('region_analysis')) if mode.get()=="特异区域 + CaSilico" else str(storage_path('runs')))
        Path(folder).mkdir(parents=True, exist_ok=True)
        os.startfile(folder)

    def close():
        if state["process"] is not None:
            messagebox.showinfo("分析正在运行", "请等待分析结束后再关闭窗口。")
        else:
            window.destroy()

    launch = ttk.Button(buttons, text="开始设计", command=start)
    launch.pack(side="left")
    ttk.Button(buttons, text="打开结果文件夹", command=open_output).pack(side="left", padx=12)
    def select_local():
        from ssu18s_import import infer_species,read_records
        files=filedialog.askopenfilenames(title='选择目标及背景18S FASTA（每个文件一个物种）',filetypes=[('FASTA','*.fa *.fasta *.fna *.fas *.fsa *.gz'),('全部','*')])
        if not files:return
        dialog=tk.Toplevel(window);dialog.title('确认物种名称，选择一个目标文件')
        canvas=tk.Canvas(dialog,width=850,height=360);canvas.pack(side='left',fill='both',expand=True)
        scroll=ttk.Scrollbar(dialog,orient='vertical',command=canvas.yview);scroll.pack(side='right',fill='y');canvas.configure(yscrollcommand=scroll.set)
        content=ttk.Frame(canvas,padding=12);canvas.create_window((0,0),window=content,anchor='nw')
        content.bind('<Configure>',lambda event:canvas.configure(scrollregion=canvas.bbox('all')))
        choice=tk.IntVar(value=0);species_vars=[]
        ttk.Label(content,text='请确认输入是核18S序列。每个文件填写一个完整物种学名，并选择目标。').grid(row=0,column=0,columnspan=3,pady=8)
        for i,file in enumerate(files):
            try:first=next(read_records(file))[0]
            except Exception as exc:messagebox.showerror('FASTA读取失败',str(exc),parent=dialog);dialog.destroy();return
            species=tk.StringVar(value=infer_species(first) or infer_species(Path(file).stem));species_vars.append(species)
            ttk.Radiobutton(content,text='目标',variable=choice,value=i).grid(row=i+1,column=0,padx=6,pady=5)
            ttk.Label(content,text=Path(file).name).grid(row=i+1,column=1,sticky='w',padx=8)
            ttk.Entry(content,textvariable=species,width=38).grid(row=i+1,column=2,padx=8)
        def confirm():
            if len(files)<2 or any(len(v.get().split())<2 for v in species_vars):
                messagebox.showerror('需要物种信息','至少选择一个目标文件和一个背景文件，并填写各自完整学名。',parent=dialog);return
            local_inputs[:]=[{'path':file,'species':v.get().strip(),'target':i==choice.get()} for i,(file,v) in enumerate(zip(files,species_vars))]
            name.set(species_vars[choice.get()].get().strip());mode.set('18S本地FASTA');coverage_min.set('');coverage_max.set('')
            append(f'已选择{len(files)}个本地文件，目标文件为{Path(files[choice.get()]).name}。\n')
            dialog.destroy()
        ttk.Button(content,text='使用这些文件',command=confirm).grid(row=len(files)+1,column=0,columnspan=3,pady=12)
    ttk.Button(buttons,text='导入FASTA并选目标文件',command=select_local).pack(side='left',padx=8)
    window.protocol("WM_DELETE_WINDOW", close)
    poll()
    window.mainloop()


if __name__ == "__main__":
    main()
