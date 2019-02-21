import urllib, urllib2

protein_names=[]
kinase_info=[]
kinases_list=[]
desired_kinases=[]
residue_info={"Phosphoserine": "S", "Phosphotyrosine": "Y", "Phosphothreonine":"T"}
# A list of expressions to remove from UniProt PTM data is created here.
# UniProt entries that do not adhere to the universal naming format are therefore omitted later on.
phrase_list=["and", "or", "alternate", "transient", "in", "vitro", "dsDNA", "kinase", "viral", "VacV", "B1",
             "SRC-type", "Tyr-kinases"]

# Remove target proteins supported by less than two peptides.
CSV = open("protein_list.csv", "r")
header = CSV.readline()
for line in CSV:
    line = line.rstrip(" \n\r")
    line = line.split(",")
    if int(line[2]) < 2:
        continue
    else: 
        protein_names.append(line[1])
CSV.close()

TSV = open("target_kinase_table.tsv", "a")
header="Target" + " " + "Position" + " " + "Residue" + " " + "Kinase" + "\n"
TSV.write(header)

# Retrieve PTM information for each target protein via UniProt API.
for entry_name in protein_names:
    target=entry_name
    kinase_info=[]
    url_req = urllib2.Request("https://www.uniprot.org/uniprot/?query=organism:9606+AND+" + entry_name + "&columns=reviewed,feature(MODIFIED%20RESIDUE)&format=tab")
    url_open = urllib2.urlopen(url_req)
    url_read = url_open.read()
    url_read = url_read.split("MOD_RES")
    for element in url_read:
        if "by" in element and "Phospho" in element:
            kinase_info.append(element)
    for entry in kinase_info:
        kinases_list=[]
        desired_kinases=[]
        entry=entry.rstrip("\n")
        entry=entry.split(" ")
        index_by=entry.index("by")
        kinases_list.extend(entry[index_by+1:])
        while "" in kinases_list:
            kinases_list.remove("")
        for x in kinases_list:
            x=x.rstrip(".,;:")
            if "ECO" not in x:
                desired_kinases.append(x)
        for phrase in phrase_list:
            while phrase in desired_kinases:
                desired_kinases.remove(phrase)
        if desired_kinases!=[]:
            for kinase in desired_kinases:
                position=entry[1]
                residue=entry[3].rstrip(';')
                if kinase=="autocatalysis":
                    auto_index=desired_kinases.index(kinase)
                    desired_kinases[auto_index]=target.rstrip("_HUMAN")
                    row=target + " " + position + " " + residue_info[residue] + " " + desired_kinases[auto_index] + "\n"
                    TSV.write(row)
                else:
                    row=target + " " + position + " " + residue_info[residue] + " " + kinase + "\n"
                    TSV.write(row)
        else:
            continue   
TSV.close()

# A list of kinases ranked by the number of phospohorylated phosphosites is created in the steps below.
all_kinases=[]
kinase_ranks=[]
top_kinases=[]

# All protein targets also extracted here for the next question.
all_targets=[]

read_TSV=open("target_kinase_table.tsv", "r")
titles=read_TSV.readline()
for line in read_TSV:
    line=line.rstrip("\n")
    line=line.split(" ")
    all_kinases.append(line[3])
    all_targets.append(line[0])
read_TSV.close()

# First extract a set of unique kinases to calculate the rankings.
unique_kinases=set(all_kinases)

for y in unique_kinases:
    counter=0
    for x in all_kinases:
        if y==x:
            counter=counter+1
    kinase_ranks.append((counter, y))

# Produce a final top kinase list according to the rank of each kinase. 
# To reflect the true count, synonymous kinases joined by a "/" are not yet separated at this stage.
kinase_ranks.sort(reverse=True)
top_kinases=[y for x,y in kinase_ranks]
print(top_kinases)

# Some kinases are corrected here.
#(1) Double kinases joined with a "/" symbol are synonymous and need separating as only one can be mapped to 
# a UniProt ID. Each time these kinases are separated, a copy of the the corresponding target is inserted 
# at an appropriate position in the targets list so that the corrected kinases map exactly to their targets.
#(2) -alpha and -beta extension in GSK3 needs converting to "A" and "B" for correct ID mapping.

corrections=[]

# Index counter determines where to insert a corresponding target when kinases have been split. 
index_counter= -1
for kinase in all_kinases:
    if "/" not in kinase and "-" not in kinase:
        corrections.append(kinase)
        index_counter+=1
    elif "/" in kinase:        
        kinase=kinase.split("/")
        corrections.extend(kinase)
        index_counter+=2
        all_targets.insert(index_counter-1, all_targets[index_counter-1])
    elif "-alpha" in kinase:
        kinase=kinase.replace("-alpha", "A")
        corrections.append(kinase)
        index_counter+=1
    elif "-beta" in kinase:
        kinase=kinase.replace('-beta', "B")
        corrections.append(kinase)
        index_counter+=1

# Map each kinase name to its corresponding *_HUMAN UniProt ID.
id_list=[]
url = 'https://www.uniprot.org/uniprot/'

for kinase in corrections:
    params = {
    'format':'tab',
    'query':'gene_exact:' + kinase + ' AND organism:homo_sapiens AND reviewed:yes',
    'columns': 'id,entry_name,genes'
    }

    data = urllib.urlencode(params)
    request = urllib2.Request(url, data)
    contact = "justynagredecka@outlook.com" 
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urllib2.urlopen(request)
    heading = response.readline()
    entries=response.read()
    new_entries=entries.split("\t")
    
    while "" in new_entries:
        new_entries.remove("")
    if new_entries!=[]:
        if "_HUMAN" in new_entries[1]:
            id_list.append(new_entries[1])
        elif "_HUMAN" not in new_entries[1]:
            id_list.append("None")
    else:
        id_list.append("None")

# Produce a unique set of interactions for Cytoscape (duplicates removed).
inter_list=[]
for n in range(0, len(id_list)):
    if id_list[n] =="None":
        continue
    else:
        inter_list.append((id_list[n], all_targets[n]))
unique_inter=set(inter_list)

# *FULL NETWORK SIF file (too coomplex to visualise graphically).
SIF_FILE=open("full_network.sif", "a")
for interaction in unique_inter:
    row = interaction[0] + " " + "phosphorylates" + " " + interaction[1] + "\n"
    SIF_FILE.write(row)
SIF_FILE.close()  

# Create a simplified SIF file for kinase proteins only.
unique_inter_list=list(unique_inter)
kinases_only=[x for x,y in unique_inter_list]
SIF_FILE=open("kinase_network.sif", "a")
for interaction in unique_inter_list:
    if interaction[1] not in kinases_only:
        continue
    else:
        row = interaction[0] + " " + "phosphorylates" + " " + interaction[1] + "\n"
        SIF_FILE.write(row)
        
SIF_FILE.close()
