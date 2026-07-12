with open("/vol8/home/xtu_pcy/l_s/PerfBench-BUAAHPC/PERFBENCH_USAGE.md", "r") as f:
    text = f.read()

target = "---\n\n## 测试示例"

replace = """---\n\n## 已修改的 PerfBench 源码\n\n### 1. `perfbench/__main__.py`\n\n在 `generate_certificate_for_test` 调用前加了 `time.sleep(3)`，等待监控脚本写完最终的 sacct 日志。\n\n### 2. `perfbench/utils/result_handler.py`\n\n**`parse_sacct` 方法**：加了 `sacct_files.sort()`，确保按文件名（时间戳）顺序处理日志文件。\n\n**`get_elapsed_time_slurm` 方法**：改为从后往前找第一条 `COMPLETED` 或 `FAILED` 状态且有 `Elapsed` 值的记录，避免最后一条记录状态为 `RUNNING` 或 `Elapsed` 为空时读到错误值。\n\n---\n\n## 测试示例"""

if target in text:
    text = text.replace(target, replace)
    with open("/vol8/home/xtu_pcy/l_s/PerfBench-BUAAHPC/PERFBENCH_USAGE.md", "w") as f:
        f.write(text)
    print("REPLACED")
else:
    print("NOT FOUND")
    idx = text.find("---\n\n## 测试示例")
    print(f"Found at {idx}")
