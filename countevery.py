import json 
  

# JSON file 
f = open ('wholecombine.json', "r") 
  
# Reading from file 
data = json.loads(f.read()) 
  
# Iterating through the json 
# list 
#print(data["Sample"])
 
wholelistmut1={}
for x in data["Sample"]:
    for y in data["Sample"][x]["VariantCalls"]:
        if y not in wholelistmut1:
            wholelistmut1[y]=1
        else: wholelistmut1[y]=wholelistmut1[y]+1

print(wholelistmut1)
# Closing file 
f.close() 