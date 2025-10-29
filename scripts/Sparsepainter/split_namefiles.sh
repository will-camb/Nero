input_file="target_samples"
output_dir="namefile"
mkdir -p "$output_dir"

lines_per_file=1000
file_count=1
line_count=0

while IFS= read -r line; do
    if [ $line_count -eq 0 ]; then
        output_file="${output_dir}/target${file_count}.txt"
        > "$output_file"
    fi

    echo "$line" >> "$output_file"

    ((line_count++))

    if [ $line_count -eq $lines_per_file ]; then
        ((file_count++))
        line_count=0
        echo $file_count
    fi
done < "$input_file"