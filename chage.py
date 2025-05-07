def process_data(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # 跳过空行
            if not line.strip():
                outfile.write(line)
                continue

            # 分割每行数据
            parts = line.split()
            if len(parts) != 7:
                outfile.write(line)  # 如果不是7列数据，保持原样
                continue

            try:
                # 提取各列数据
                x = parts[0]
                y = parts[1]
                z = parts[2]
                sp = float(parts[3])
                spx = float(parts[4])
                spy = float(parts[5])

                # 计算新值
                new_spx = sp * spx
                new_spy = sp * spy

                # 写入新格式的行
                outfile.write(f"{x} {y} {z} {new_spx:.15g} {new_spy:.15g}\n")
            except ValueError:
                outfile.write(line)  # 如果转换失败，保持原样


if __name__ == "__main__":
    input_file = "4-28foce.txt"
    output_file = "4-28foce_processed.txt"

    process_data(input_file, output_file)
    print(f"处理完成，结果已保存到 {output_file}")