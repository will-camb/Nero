from ms_output_to_cp_input import RunMsWithCpOutput

model = RunMsWithCpOutput(nhaps=[2,2,2,2,2], number_recipient_haps=2)

#model = RunMsWithFsOutput(nhaps=[8, 4, 4, 4, 4], length=1982950)
model.run()
