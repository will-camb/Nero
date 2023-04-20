# https://stackoverflow.com/questions/26124417/how-to-convert-a-csv-file-to-parquet
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-chr",
                    help="Chromosome number",
                    required=True)
args = parser.parse_args()

for ancestry in ['African', 'Yamnaya', 'Farmer', 'CHG', 'WHG', 'EHG', 'EastAsian']:
    csv_file = f'{ancestry}.{args.chr}.master_all_copyprobsperlocus.txt.gz'
    parquet_file = f'parquet_files/{ancestry}.{args.chr}.master_all_copyprobsperlocus.test.parquet'
    chunksize = 1000

    csv_stream = pd.read_csv(csv_file, chunksize=chunksize, low_memory=False, sep=" ")

    for i, chunk in enumerate(csv_stream):
        print("Chunk", i)
        if i == 0:
            # Guess the schema of the CSV file from the first chunk
            parquet_schema = pa.Table.from_pandas(df=chunk).schema
            # Open a Parquet file for writing
            parquet_writer = pq.ParquetWriter(parquet_file, parquet_schema, compression='snappy')
        # Write CSV chunk to the parquet file
        table = pa.Table.from_pandas(chunk, schema=parquet_schema)
        parquet_writer.write_table(table)

parquet_writer.close()
