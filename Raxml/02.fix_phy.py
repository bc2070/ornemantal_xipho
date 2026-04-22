import sys
import os

def main():
    if len(sys.argv) < 2:
        print("Usage: python fix_phy.py input.phy")
        sys.exit(1)

    input_file = sys.argv[1]
    # 生成输出文件名，例如 input_fixed.phy
    base, ext = os.path.splitext(input_file)
    output_file = f"{base}_fixed{ext}"

    sequences = []
    
    print(f"Reading {input_file}...")

    with open(input_file, 'r') as f:
        # 读取第一行（旧的 Header），但我们不使用其中的数字，直接跳过
        first_line = f.readline()
        
        # 逐行读取剩下的序列数据
        for line_num, line in enumerate(f, start=2):
            line = line.strip()
            if not line:
                continue
            
            # 使用 split(maxsplit=1) 分割名字和序列
            # 这样可以处理名字和序列之间有多个空格或制表符的情况
            parts = line.split(maxsplit=1)
            
            if len(parts) < 2:
                print(f"Warning: Line {line_num} seems malformed (no sequence found). Skipping.")
                continue
            
            name = parts[0]
            seq = parts[1]
            
            # 确保序列中没有空格或制表符（防止格式混乱）
            seq = seq.replace(" ", "").replace("\t", "")
            
            sequences.append((name, seq))

    if not sequences:
        print("Error: No sequences found in the file.")
        sys.exit(1)

    # 获取序列数量
    num_seqs = len(sequences)
    
    # 获取第一条序列的长度作为基准
    ref_len = len(sequences[0][1])
    
    print(f"Detected: {num_seqs} sequences.")
    print(f"Base length (from first seq): {ref_len}")

    # 检查所有序列长度是否一致（PHY文件必须要求等长）
    # 如果你删除了 * 导致长度不一致，这里会报错提醒
    for name, seq in sequences:
        if len(seq) != ref_len:
            print(f"\n[ERROR] Sequence length mismatch detected!")
            print(f"Sample '{sequences[0][0]}' length: {ref_len}")
            print(f"Sample '{name}' length: {len(seq)}")
            print("Phylip (.phy) files required all sequences to be of equal length.")
            print("Please check your data processing step (e.g. removal of *).")
            sys.exit(1)

    # 写入新文件
    with open(output_file, 'w') as f_out:
        # 1. 写入修正后的 Header
        f_out.write(f"{num_seqs} {ref_len}\n")
        
        # 2. 写入序列
        for name, seq in sequences:
            f_out.write(f"{name}\t{seq}\n")

    print(f"\nSuccess! Fixed file saved to: {output_file}")
    print(f"Header updated to: {num_seqs} {ref_len}")

if __name__ == "__main__":
    main()
