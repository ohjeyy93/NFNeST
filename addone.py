with open("genesprefix.gff", "r") as f1:
    with open("genesfix.gff", "w") as f2:
        for line in f1:
            if line.startswith("##") == False and line.find("extracted")==-1 and line.find("region")==-1:
                templine=""
                count=0
                #print(line.count("\t"))
                for word in line.split():
                    if count!=3: templine+=word+"\t"
                    if count==3 and word!="region":
                        #print(word)
                        templine+=str(int(word)+1)+"\t"
                    count+=1
                templine+="\n"
                f2.write(templine)
            elif line.startswith("##") == False and line.find("extracted")!=-1 and line.find("region")!=-1:
                #print("true")
                templine=""
                count=0
                #print(line.count("\t"))
                for word in line.split():
                    if count<10:
                        templine+=word+"\t"
                    if count==4:
                        #print(word)
                        templine+=str(int(word)+1)+"\t"
                    if count==10:
                        templine+=word+"\n"
                    count+=1
                f2.write(templine)
            elif line.startswith("##"): f2.write(line)
