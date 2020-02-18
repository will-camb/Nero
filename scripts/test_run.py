from ms_output_to_fs_input import RunMsWithFsOutput

model = RunMsWithFsOutput(nhaps=[8, 4, 4, 4, 4], length=1982950)
model.run()
