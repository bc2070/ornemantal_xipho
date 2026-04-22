import sys
import numpy as np
from scipy import stats

if len(sys.argv) != 2:
    print(f"用法: python {sys.argv[0]} <input_file>")
    sys.exit(1)

input_file = sys.argv[1]
fifth_column_data = []

try:
    with open(input_file, 'r') as f:
        for line in f:
            try:
                parts = line.strip().split()
                if len(parts) >= 5:
                    value = float(parts[4])
                    fifth_column_data.append(value)
            except (ValueError, IndexError):
                continue
except FileNotFoundError:
    print(f"错误: 文件 '{input_file}' 未找到。")
    sys.exit(1)

if len(fifth_column_data) < 2:
    print("错误: 数据不足，无法进行分析。")
    sys.exit(1)

data_array = np.array(fifth_column_data)

percentile_2_5 = np.percentile(data_array, 2.5)
percentile_97_5 = np.percentile(data_array, 97.5)

print(f"2.5分位数: {percentile_2_5}")
print(f"97.5分位数: {percentile_97_5}")
print("\n检验数据分布峰值与0的差异 (Bootstrap CI for Mode):")

def estimate_mode(data):
    if np.std(data) == 0:
        return data[0]
    kde = stats.gaussian_kde(data)
    x_grid = np.linspace(min(data), max(data), 1000)
    density = kde.evaluate(x_grid)
    return x_grid[np.argmax(density)]

original_mode = estimate_mode(data_array)
print(f"数据分布的估计峰值 (Mode): {original_mode}")

n_bootstrap = 10000
bootstrap_modes = []
for _ in range(n_bootstrap):
    sample = np.random.choice(data_array, size=len(data_array), replace=True)
    try:
        mode = estimate_mode(sample)
        bootstrap_modes.append(mode)
    except np.linalg.LinAlgError:
        continue

if not bootstrap_modes:
    print("无法通过自助法计算置信区间。")
    sys.exit(1)

ci_lower = np.percentile(bootstrap_modes, 2.5)
ci_upper = np.percentile(bootstrap_modes, 97.5)

print(f"峰值的95%置信区间: [{ci_lower}, {ci_upper}]")

if ci_lower <= 0 <= ci_upper:
    print("结论: 峰值与0无显著差异 (0在95%置信区间内)。")
else:
    print("结论: 峰值与0有显著差异 (0不在95%置信区间内)。")
