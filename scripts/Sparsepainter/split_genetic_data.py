import gzip

chrom = '6'
input_file = f"target.{chrom}.phase"
lines_per_file = 2000


def split_file(input_file_, lines_per_file_, chr_):
    with open(input_file_, 'r') as file:
        next(file)  # Skip the first line
        second_line = next(file)
        third_line = next(file)
    with open(input_file_, 'r') as file:
        for _ in range(3):
            next(file)
        file_number = 1
        line_count = 0
        output_file = gzip.open(f"chr{chr_}split/chr{chr_}_target{file_number}.phase", 'wt')
        output_file.write(f"{lines_per_file_}\n")
        output_file.write(second_line)
        output_file.write(third_line)
        for line in file:
            if line_count == lines_per_file_:
                output_file.close()
                file_number += 1
                line_count = 0
                output_file = gzip.open(f"chr{chr_}split/chr{chr_}_target{file_number}.phase", 'wt')
                output_file.write(f"{lines_per_file_}\n")
                output_file.write(second_line)
                output_file.write(third_line)
            output_file.write(line)
            line_count += 1
        output_file.close()


split_file(input_file, lines_per_file, chrom)
