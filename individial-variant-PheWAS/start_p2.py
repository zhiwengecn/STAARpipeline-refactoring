#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import os
import sys
import argparse
import subprocess
from datetime import datetime
from typing import List, Dict, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
import yaml


# ===================== 配置加载函数 =====================
def load_config(config_path: str) -> Dict:
    """从 YAML 文件加载配置"""
    config_file = Path(config_path)
    if not config_file.exists():
        print(f"[错误] 配置文件不存在：{config_file}")
        sys.exit(1)
    
    try:
        with open(config_file, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        return config
    except yaml.YAMLError as e:
        print(f"[错误] YAML格式错误：{e}")
        sys.exit(1)
    except Exception as e:
        print(f"[错误] 读取配置文件失败：{e}")
        sys.exit(1)


def validate_config(config: Dict) -> bool:
    """验证配置完整性"""
    required_sections = ['paths', 'tasks', 'execution']
    for section in required_sections:
        if section not in config:
            print(f"[错误] 配置缺少必需的部分：{section}")
            return False
    
    # 验证路径配置
    required_paths = ['r_script', 'jobs_path', 'pheno_path', 'step1_path', 
                     'out_path', 'step2_1_path', 'step2_2_path']
    for path_key in required_paths:
        if path_key not in config['paths']:
            print(f"[错误] 配置缺少必需的路径：{path_key}")
            return False
    
    # 验证任务配置
    if not config['tasks']:
        print(f"[错误] 配置中至少需要定义一个任务")
        return False
    
    for i, task in enumerate(config['tasks']):
        required_fields = ['begin', 'end', 'stage', 'arrayid']
        for field in required_fields:
            if field not in task:
                print(f"[错误] 任务 {i+1} 缺少必需字段：{field}")
                return False
    
    return True


def read_task_configs(config_data: Dict) -> List[Dict[str, int]]:
    """从配置数据中提取任务配置"""
    tasks = []
    for task_config in config_data.get('tasks', []):
        task_params = {
            'begin': task_config['begin'],
            'end': task_config['end'],
            'stage': task_config['stage'],
            'arrayid': task_config['arrayid']
        }
        tasks.append(task_params)
    return tasks


def read_subset_num(subset_num_file: Path) -> Optional[int]:
    """从 .txt 文件读取 subset_num（假设文件内容为单个整数）"""
    if not subset_num_file.exists():
        print(f"[错误] subset_num 文件不存在：{subset_num_file}")
        return None
    try:
        content = subset_num_file.read_text(encoding="utf-8").strip()
        num = int(content)
        if num <= 0:
            print(f"[错误] subset_num 必须为正整数，实际为 {num}，文件：{subset_num_file}")
            return None
        return num
    except Exception as e:
        print(f"[错误] 读取 subset_num 失败：{subset_num_file}，错误：{e}")
        return None


def generate_all_subtasks(tasks_config: List[Dict], step2_1_path: str) -> List[Tuple[int, int, str, int, int]]:
    """为每个任务生成 kk=1..subset_num 的子任务"""
    all_subtasks = []
    
    for task in tasks_config:
        begin = task['begin']
        end = task['end']
        stage = task['stage']
        arrayid = task['arrayid']
        
        # 构造 subset_num 文件路径
        filename = f"subsetnum_{begin}_to_{end}_array{arrayid}_{stage}.txt"
        subset_num_file = Path(step2_1_path) / filename

        subset_num = read_subset_num(subset_num_file)
        if subset_num is None:
            print(f"[跳过] 无法获取 subset_num，跳过任务组：{begin}-{end} {stage} {arrayid}")
            continue

        print(f"[INFO] {begin}-{end} {stage} {arrayid} -> subset_num = {subset_num}, 将生成 {subset_num} 个 kk 子任务")

        for kk in range(1, subset_num + 1):
            all_subtasks.append((begin, end, stage, arrayid, kk))

    return all_subtasks


def build_r_command(begin: int, end: int, stage: str, arrayid: int, kk: int, 
                   paths: Dict) -> List[str]:
    """构造 R 命令行调用"""
    return [
        "R",
        "--slave",
        "--no-restore",
        "--file=" + paths['r_script'],
        "--args",
        "--jobs_path", paths['jobs_path'],
        "--log_path", "",  # 运行时动态设置
        "--pheno_path", paths['pheno_path'],
        "--step1_path", paths['step1_path'],
        "--out_path", paths['out_path'],
        "--step2.1_path", paths['step2_1_path'],
        "--step2.2_path", paths['step2_2_path'],
        "--arrayid", str(arrayid),
        "--begin", str(begin),
        "--end", str(end),
        "--stage", stage,
        "--kk", str(kk)
    ]


def run_single_task(task: Tuple[int, int, str, int, int], log_root: Path, 
                   paths: Dict, dry_run: bool) -> dict:
    """运行单个 R 任务（包含 kk），记录日志与状态"""
    begin, end, stage, arrayid, kk = task
    task_id = f"{begin}-{end}_{stage}_{arrayid}_kk{kk}"
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    log_filename = f"r_task_{task_id}_{timestamp}.log"
    log_file = log_root / log_filename

    cmd = build_r_command(begin, end, stage, arrayid, kk, paths)
    # 动态设置日志路径参数
    cmd[7] = "--log_path"
    cmd[8] = str(log_root)

    result = {
        "task": task,
        "task_id": task_id,
        "cmd": " ".join(cmd),
        "log_file": str(log_file),
        "start_time": datetime.now(),
        "end_time": None,
        "returncode": None,
        "success": False,
        "error": None
    }

    if dry_run:
        print(f"[DRY-RUN] Would run: {result['cmd']} > {log_file}")
        result["success"] = True
        result["end_time"] = datetime.now()
        return result

    try:
        with log_file.open("w", encoding="utf-8") as lf:
            lf.write(f"Command: {' '.join(cmd)}\n")
            lf.write(f"Task ID: {task_id}\n")
            lf.write(f"Parameters: begin={begin}, end={end}, stage={stage}, arrayid={arrayid}, kk={kk}\n")
            lf.write(f"Start: {result['start_time']}\n\n")
            lf.flush()

            proc = subprocess.Popen(
                cmd,
                stdout=lf,
                stderr=subprocess.STDOUT,
                cwd=None,
                env=os.environ
            )
            retcode = proc.wait()
            result["returncode"] = retcode
            result["success"] = (retcode == 0)
    except Exception as e:
        result["error"] = str(e)
    finally:
        result["end_time"] = datetime.now()

    return result


def main():
    parser = argparse.ArgumentParser(description="并发提交多个 R 任务（支持 kk 子任务）")
    parser.add_argument("--config", default="config.yaml", help="YAML 配置文件路径")
    parser.add_argument("--max-jobs", type=int, help="最大并发任务数（覆盖配置）")
    parser.add_argument("--dry-run", action="store_true", help="仅打印命令，不实际运行")
    parser.add_argument("--validate-only", action="store_true", help="仅验证配置文件，不执行任务")
    args = parser.parse_args()

    # Step 1: 加载配置
    print(f"📥 从配置文件加载：{args.config}")
    config = load_config(args.config)
    
    # 验证配置
    if not validate_config(config):
        print("❌ 配置验证失败")
        sys.exit(1)

    # 合并命令行参数
    if args.max_jobs:
        config['execution']['max_concurrent_jobs'] = args.max_jobs
    if args.dry_run:
        config['execution']['dry_run'] = True

    print(f"📋 配置加载成功")
    print(f"📁 日志根目录: {config['paths']['log_root']}")
    print(f"🔄 最大并发数: {config['execution']['max_concurrent_jobs']}")
    print(f"🎯 任务数量: {len(config['tasks'])}")

    # Step 2: 提取任务配置
    config_tasks = read_task_configs(config)
    print(f"📥 共读取到 {len(config_tasks)} 个主任务")

    # 如果只是验证模式，到此结束
    if args.validate_only:
        print("✅ 配置验证完成")
        sys.exit(0)

    # Step 3: 展开为所有 kk 子任务
    all_subtasks = generate_all_subtasks(config_tasks, config['paths']['step2_1_path'])
    if not all_subtasks:
        print("❌ 没有生成任何有效的子任务（可能 subset_num 文件缺失或无效）")
        sys.exit(1)

    print(f"📦 共展开为 {len(all_subtasks)} 个 R 子任务（含 kk）")

    # Step 4: 创建日志目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_root = Path(config['paths']['log_root']) / f"r_run_{timestamp}"
    log_root.mkdir(parents=True, exist_ok=True)
    print(f"📂 日志将保存至：{log_root}")

    # Step 5: 并发执行
    success = 0
    failed = 0
    total_tasks = len(all_subtasks)
    
    max_jobs = config['execution']['max_concurrent_jobs']
    dry_run = config['execution']['dry_run']
    
    print(f"🚀 开始执行 {total_tasks} 个任务，最大并发数：{max_jobs}")
    print("=" * 80)

    start_time = datetime.now()
    
    if dry_run:
        print("\n🧪 干跑模式 - 以下命令将被执行:")
        for task in all_subtasks:
            begin, end, stage, arrayid, kk = task
            cmd = build_r_command(begin, end, stage, arrayid, kk, config['paths'])
            cmd[7] = "--log_path"
            cmd[8] = str(log_root)
            print(f"CMD: {cmd[0]} {' '.join(cmd[1:])}")
        return

    with ProcessPoolExecutor(max_workers=max_jobs) as executor:
        futures = {
            executor.submit(run_single_task, task, log_root, config['paths'], dry_run): task
            for task in all_subtasks
        }

        for future in as_completed(futures):
            res = future.result()
            status = "✅ 成功" if res["success"] else "❌ 失败"
            if res["error"]:
                status += f" (异常: {res['error']})"
            
            duration = res["end_time"] - res["start_time"]
            
            print(f"{status} | {res['task_id']} | 耗时: {duration} | 日志: {Path(res['log_file']).name}")

            if res["success"]:
                success += 1
            else:
                failed += 1

    end_time = datetime.now()
    total_duration = end_time - start_time

    print("=" * 80)
    print(f"📊 执行总结：")
    print(f"   成功: {success}")
    print(f"   失败: {failed}")
    print(f"   总计: {total_tasks}")
    print(f"   总耗时: {total_duration}")
    
    if failed > 0:
        print("⚠️  有任务失败，请检查对应日志文件排查问题")
        sys.exit(1)
    else:
        print("🎉 所有任务执行成功！")


if __name__ == "__main__":
    main()
