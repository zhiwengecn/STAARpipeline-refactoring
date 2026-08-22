#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MultiSTAAR Pipeline Launcher with Concurrent kk Tasks
Submit multiple STAARpipeline R tasks concurrently with kk sub-task ranges
"""

from pathlib import Path
import os
import sys
import yaml
import argparse
import subprocess
from datetime import datetime
from typing import List, Dict, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed


# ===================== Configuration Loader =====================
def load_config(config_path: Path) -> Dict:
    """Load configuration from YAML file"""
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        
        meta_summary_path = config.get('paths', {}).get('meta_summary')
        if meta_summary_path and Path(meta_summary_path).exists():
            subset_variants = {}
            with open(meta_summary_path, 'r', encoding='utf-8') as msf:
                header = next(msf).strip().split('\t')
                subset_index = header.index("subset_variants_num")
                for line in msf:
                    parts = line.strip().split('\t')
                    array_num = int(parts[0])
                    subset_variants[array_num] = int(parts[subset_index])
            config['subset_variants'] = subset_variants
        else:
            print(f"[WARNING] Analysis_Meta_Summary.txt not found or not configured: {meta_summary_path}")
            config['subset_variants'] = {}
        return config
    except Exception as e:
        print(f"[ERROR] Failed to load config file: {e}")
        sys.exit(1)

# =====================================================================

def read_task_configs(config_data: Dict) -> List[Dict[str, int]]:
    """Extract task configurations from config data"""
    tasks = []
    for task_config in config_data.get('tasks', []):
        task_params = {
            'array_id': task_config['array_id'],
            'save_chunk_size': task_config['save_chunk_size'],
            'kk_start': task_config['kk_start'],
            'kk_end': task_config['kk_end']
        }
        tasks.append(task_params)
    return tasks


def build_r_command(params: Dict[str, int], config_paths: Dict, config: Dict) -> List[str]:
    array_id = params['array_id']
    subset_num = config['subset_variants'].get(array_id)

    return [
        "R",
        "--slave",
        "--no-restore",
        "--file=" + config_paths['r_script'],
        "--args",
        "--adgs", config_paths['adgs'],
        "--null", config_paths['null_model'],
        "--jobs", config_paths['jobs'],
        "--out", config_paths['out_prefix'],
        "--p1_path", config_paths['p1_path'],
        "--p2_path", config_paths['p2_path'],
        "--log_path", config_paths['log_root'],
	    "--array_id", str(params['array_id']),
        "--save_chunk", str(params['save_chunk_size']),
        "--subset_num", str(subset_num),
        "--kk_start", str(params['kk_start']),
        "--kk_end", str(params['kk_end'])
    ]


def run_single_task(params: Dict[str, int], log_root: Path, dry_run: bool, config_paths: Dict, config: Dict) -> dict:
    """Run a single R task and record logs and status"""
    task_id = f"arr{params['array_id']}_kk{params['kk_start']}-{params['kk_end']}"
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    log_filename = f"staar_task_{task_id}_{timestamp}.log"
    log_file = log_root / log_filename

    cmd = build_r_command(params, config_paths, config)

    result = {
        "params": params,
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
            lf.write(f"Parameters: {params}\n")
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


def validate_parameters(params: Dict[str, int]) -> bool:
    """Validate task parameters"""
    # 检查 array 范围
    if params['array_id'] < 1 or params['array_id'] > 288:
        print(f"[ERROR] array_id must be within 1–288, current: {params['array_id']}")
        return False
    
    # 检查 kk 范围
    if params['kk_start'] < 1 or params['kk_end'] < 1:
        print(f"[ERROR] kk_start and kk_end must be greater than 0")
        return False
    if params['kk_end'] < params['kk_start']:
        print(f"[ERROR] kk_end ({params['kk_end']}) must be >= kk_start ({params['kk_start']})")
        return False
    
    # 检查 save_chunk_size
    if params['save_chunk_size'] < 1:
        print(f"[ERROR] save_chunk_size must be greater than 0")
        return False
    
    return True


def main():
    parser = argparse.ArgumentParser(description="Submit multiple STAARpipeline R tasks concurrently")
    parser.add_argument("--config", default="config.yaml", help="Path to YAML config file")
    parser.add_argument("--max-jobs", type=int, help="Maximum number of concurrent jobs (overrides config)")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    parser.add_argument("--validate-only", action="store_true", help="Validate config only, do not submit jobs")
    args = parser.parse_args()

    # Step 1: 加载配置
    config_path = Path(args.config)
    if config_path.exists():
        print(f"📥 Loading configuration from: {config_path}")
        config = load_config(config_path)
    else:
        print(f"⚠️ Config file not found, using default path: {config_path}")
        sys.exit(1)

    # 合并命令行参数
    if args.max_jobs:
        config['execution']['max_concurrent_jobs'] = args.max_jobs
    if args.dry_run:
        config['execution']['dry_run'] = True

    # Step 2: 提取任务配置
    config_tasks = read_task_configs(config)
    if not config_tasks:
        print("❌ No valid task configurations found. Please check your config file.")
        sys.exit(1)

    print(f"📥 Loaded {len(config_tasks)} task configurations")

    # Step 3: 验证参数
    valid_tasks = []
    invalid_count = 0
    for i, params in enumerate(config_tasks, 1):
        if validate_parameters(params):
            valid_tasks.append(params)
        else:
            invalid_count += 1
            print(f"[ERROR] Task #{i} has invalid parameters and will be skipped")

    if invalid_count > 0:
        print(f"⚠️  Found {invalid_count} invalid task configuration(s)")

    if not valid_tasks:
        print("❌ No valid tasks to execute")
        sys.exit(1)

    print(f"✅ Valid tasks: {len(valid_tasks)}")

    # 如果只是验证模式，到此结束
    if args.validate_only:
        print("✅ Parameter validation completed")
        sys.exit(0)

    # Step 4: 创建日志目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_root = Path(config['paths']['log_root']) / f"staar_run_{timestamp}"
    log_root.mkdir(parents=True, exist_ok=True)
    print(f"📂 Logs will be saved to: {log_root}")

    # Step 5: 并发执行
    success = 0
    failed = 0
    total_tasks = len(valid_tasks)
    
    max_jobs = config['execution']['max_concurrent_jobs']
    dry_run = config['execution']['dry_run']
    
    print(f"🚀 Submitting {total_tasks} tasks with max concurrency: {max_jobs}")
    print("=" * 80)

    start_time = datetime.now()
    
    with ProcessPoolExecutor(max_workers=max_jobs) as executor:
        futures = {
            executor.submit(run_single_task, task, log_root, dry_run, config['paths'], config): task
            for task in valid_tasks
        }

        for future in as_completed(futures):
            res = future.result()
            status = "✅ Success" if res["success"] else "❌ Failed"
            if res["error"]:
                status += f" (Exception: {res['error']})"
            
            duration = res["end_time"] - res["start_time"]
            
            print(f"{status} | {res['task_id']} | Duration: {duration} | Log: {Path(res['log_file']).name}")

            if res["success"]:
                success += 1
            else:
                failed += 1

    end_time = datetime.now()
    total_duration = end_time - start_time

    print("=" * 80)
    print(f"📊 Execution Summary:")
    print(f"   Success: {success}")
    print(f"   Failed: {failed}")
    print(f"   Total: {total_tasks}")
    print(f"   Runtime: {total_duration}")
    
    if failed > 0:
        print("⚠️  Some tasks failed. Please check the corresponding log files.")
        sys.exit(1)
    else:
        print("🎉 All tasks completed successfully!")


if __name__ == "__main__":
    main()
