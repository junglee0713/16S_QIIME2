def get_sample(sample_fp):
   with open(sample_fp) as f:
      lines = f.readlines()
   samples = []
   for line in lines:
      samples.append(line.split("\t")[0])
   if "SampleID" in samples: samples.remove("SampleID")
   if "#SampleID" in samples: samples.remove("#SampleID") 
   return(samples)
