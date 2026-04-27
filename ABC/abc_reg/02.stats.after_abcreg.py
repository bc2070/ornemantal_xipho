import sys
import numpy as np
from scipy import stats

if len(sys.argv) != 2:
    print(f"python {sys.argv[0]} <input_file>")
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
    print(f"Error: file '{input_file}' not found。")
    sys.exit(1)

if len(fifth_column_data) < 2:
    print("Error: Insufficient data, unable to conduct analysis.")
    sys.exit(1)

data_array = np.array(fifth_column_data)

percentile_2_5 = np.percentile(data_array, 2.5)
percentile_97_5 = np.percentile(data_array, 97.5)

print(f"percentile_2.5: {percentile_2_5}")
print(f"percentile_97.5: {percentile_97_5}")
print("\n(Bootstrap CI for Mode):")

def estimate_mode(data):
    if np.std(data) == 0:
        return data[0]
    kde = stats.gaussian_kde(data)
    x_grid = np.linspace(min(data), max(data), 1000)
    density = kde.evaluate(x_grid)
    return x_grid[np.argmax(density)]

original_mode = estimate_mode(data_array)
print(f"(Mode): {original_mode}")

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
    print("The confidence interval cannot be calculated using the self-service method.")
    sys.exit(1)

ci_lower = np.percentile(bootstrap_modes, 2.5)
ci_upper = np.percentile(bootstrap_modes, 97.5)

print(f"95% confidence interval for the peak value: [{ci_lower}, {ci_upper}]")

if ci_lower <= 0 <= ci_upper:
    print("no significant difference")
else:
    print("significant difference")
