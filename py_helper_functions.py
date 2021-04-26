import os

def parseMegan(filename, prefix=""):
    ''' Takes the MEGAN_info file generated from the MEGAN GUI and split it into the
        respective categories (TAX, INTERPRO2GO etc). '''
    output = {}
    key = ""
    data = ""
    with open(filename,"r") as f:
        while True:
            line = f.readline().strip()
            if line == "END_OF_DATA_TABLE":
                break
            elif line.split("\t")[0] == "@Names":
                data = line[6:]
                data = "CATEGORY\tNUM" + data
            elif line[0] == "@":
                continue
            else:
                key = line.split("\t")[0]
                if key not in output.keys():
                    output[key] = []
                output[key].append(line)
    for key, value in output.items():
        file = prefix + "_" + key + ".tsv"
        with open(file, "w") as newfile:
            newfile.write(data+"\n")
            newfile.write('\n'.join(line for line in value))

def interproscan_reformat(filename):
    ''' Reformat the INTERPROSCAN to GO mapping to be more consistent and easier for downstream analysis.'''
    interpro_id_ls =[]
    interpro_name_ls = []
    go_id_ls = []
    go_name_ls = []
    interpro, go = "", ""
    
    with open(filename, "r") as f:
        for line in f.readlines():
            if line[0] == "!":
                continue
            else:
                interpro, go = line.split(" > ")
                # interpro processing
                interpro = interpro.split()
                interpro_id = interpro[0].split(":")[1]
                interpro_id_ls.append(interpro_id.strip())
                interpro_name = " ".join(interpro[1:])
                interpro_name_ls.append(interpro_name.strip())
                # go processing
                go_name, go_id = go.split(" ; ")
                go_id_ls.append(go_id.strip())
                go_name_ls.append(go_name.strip())

    newfile_name = "INTERPRO2GO_MAP_CLEANED.tsv"
    with open(newfile_name, "w") as newfile:
        for a,b,c,d in zip(interpro_id_ls,interpro_name_ls,go_id_ls,go_name_ls):
            newfile.write("\t".join([a,b,c,d])+"\n")

def interproscan_goatools(filename, output="interproscan_goatools.txt"):
    mapping_data = {}
    with open(filename, "r") as f:
        for line in f.readlines():
            line = line.split("\t")
            if line[0][3:] not in mapping_data.keys():
                mapping_data[line[0][3:]] = []
            mapping_data[line[0][3:]].append(line[2])
    with open(output, "w") as out:
        for key, value in mapping_data.items():
            out.write(key+"\t"+";".join(value)+"\n")

def combine_bracken_output(filepath,level="P"):
    file_list = os.listdir(filepath)
    main_dic = {}
    #read in all data
    for file in file_list:
        sample = file.replace("_bracken_phylums.kreport", "")
        main_dic[sample] = {}
        with open(os.path.join(filepath,file), "r") as f:
            lines = f.readlines()
            for line in lines:
                line = line.split("\t")
                if line[3] == level:
                    main_dic[sample][line[5].strip()] = line[1]
    all_taxa = set()
    #get unique taxas
    for key in main_dic.keys():
        all_taxa.update(list(main_dic[key].keys()))
    out = ["taxa"]
    out.extend(all_taxa)
    for key in main_dic.keys():
        out[0] += "\t" + key
        for i in range(1, len(out)):
            taxa = out[i].split("\t")[0]
            out[i] += "\t" + main_dic[key].get(taxa, "0")
    with open("bracken_combined.tsv", "w") as f:
        f.write("\n".join(out))



#interproscan_reformat("INTERPRO2GO_MAP.txt")
#parseMegan("daa2rma.megan", "rma")
#parseMegan("root_4m_info", prefix="root4m")
#parseMegan("bulk_4m_info", prefix="bulk4m")
#interproscan_goatools("INTERPRO2GO_MAP_CLEANED.tsv")
#combine_bracken_output("C:\\Users\\YZ\\Desktop\\FYP\\dip_metagenome\\results\\bracken_kreport",level="P")
