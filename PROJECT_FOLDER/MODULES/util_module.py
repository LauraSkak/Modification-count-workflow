def create_sample_dict(sample_file):
    
    sample_dict = {}

    with open(sample_file, 'r') as file:
        lines = file.readlines()

        for line in lines[1:]:
            row = line.strip().split("\t")

            somaseqid = row[0]
            sampleid = row[1]
            
            values = [sampleid]
            
            if len(row)> 2:
                for i in range(2,len(row)):
                    samplegroup = row[i]
                    
                    values.append(samplegroup)
            
            sample_dict[somaseqid]=values
            
    return sample_dict