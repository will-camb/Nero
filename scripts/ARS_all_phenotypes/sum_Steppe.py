import pandas as pd

for i in range(1, 23):
    try:
        Yamnaya_cols = pd.read_csv("temp.Yamnaya." + str(i) + ".master_all_copyprobsperlocus.txt.gz", nrows=0).columns
        types_dict = {'ID': str}
        types_dict.update({col: 'float' for col in Yamnaya_cols if col not in types_dict})
        Yamnaya = pd.read_csv("temp.Yamnaya." + str(i) + ".master_all_copyprobsperlocus.txt.gz", dtype=types_dict)
        Yamnaya.set_index("ID", inplace=True)
        EastAsian = pd.read_csv("temp.EastAsian." + str(i) + ".master_all_copyprobsperlocus.txt.gz", dtype=types_dict)
        EastAsian.set_index("ID", inplace=True)
        Steppe = Yamnaya + EastAsian
        Steppe.to_csv("temp.Steppe." + str(i) + ".master_all_copyprobsperlocus.txt.gz")
        print("Done temp.Steppe." + str(i) + ".master_all_copyprobsperlocus.txt.gz")
    except FileNotFoundError:
        print("No files found for chr" + str(i) + ". If this is expected, ignore this message")
        continue
