import requests as req
import xmltodict as xd
import pandas as pd
import time

#input folder is where the output will be automatically saved
inputPath = str(input("Enter the path to the folder\t:\t"))
inputFile = str(input("Enter the full name of the SORTED input file\t:\t"))
df = pd.read_csv(f"{inputPath}/{inputFile}")
BasePairs = df["BP"]
Chrs = df["CHR"]
SNPs = df["SNP"]
# BasePairs and Chrs series above contain the BP positions and Chromosome numbers of ALL the SNPs in the input file.


# in case the file in question requires that you search for all the rs IDs, then the following lines provide that option.
# Else enter the number of SNPs you wanna search for.
HowManyTopSNPIDs = 0
if(str(input("Enter Y if you want to search rsIDs for all SNPs in the file.\nElse press anything else\n\t\tEnter\t\t:\t"))) == "Y":
    HowManyTopSNPIDs = len(BasePairs)
else:
    HowManyTopSNPIDs = int(input("Enter the number of top SNPs, whose SNP ID you want to fetch\t:\t"))
# options over


# slicing the above two series for the SNPs in which we are interested.
TopSNP_BasePairs = BasePairs[0:HowManyTopSNPIDs]
TopSNP_CHRs = Chrs[0:HowManyTopSNPIDs]
TopSNP_Ids = SNPs[0:HowManyTopSNPIDs]
# slicing over





# the below method accesses a webpage and extracts it's source.
# it uses the hg19 build...referred to as the "Previous build". You'll see this word in the link as well
def get_source(Link : str) -> str:
    # setting up the proxy parameters for the all powerful CGNAT. Praise be to CGNAT.
    # http_proxy = "http://pdoibale:Dec%232019@10.100.24.103:8080"
    # proxies = {"http" : http_proxy, "https" : http_proxy}
    # praises over
    TextContent = ""
    SourceContent = req.Response
    try:
        SourceContent = req.get(Link, "html.parser")
        #SourceContent = req.get(Link, "html.parser", proxies=proxies)
    except Exception as e:
        print(e)
        print("\n\n\t\tERROR as ABOVE\n\n")
        time.sleep(10)
        print("\n\n")
    TextContent = str(SourceContent.text)
    return TextContent


# below method takes up the source of an HTML webpage as queried by the loop farther below.
# it extracts the rs ID and returns it
# THIS METHOD WILL FAIL THE INSTANT NCBI decides to change certain key HTML elements.
# Get the API key to alleviate all this nonsense
def get_rsid(Source : str) -> str:
    Split1, Split2 = "",""
    SNPID = ""
    try:
        Split1 = Source.split("href=\"/snp/rs", 1)[1]
        Split2 = Split1.split("\">", 1)[0]
        SNPID = "rs"+Split2
        return SNPID
    except Exception as e:
        print("Trouble finding rsID for the above.\tAdding \"nameless\"")
        print("\n\n")
        print(e)
        print("\n\n\t\tERROR as ABOVE\n\n")
        time.sleep(10)
        SNPID = "nameless_"+str(BP)
        return SNPID



# setting up some of the variables, that we shall require outside the scope of the loop
SNPIDList = list()
NumberOfSNPS = 0
# below loop is where the main action is happening
while True:
    BP = list(TopSNP_BasePairs)[NumberOfSNPS]
    Chr = list(TopSNP_CHRs)[NumberOfSNPS]
    SNP_Ids = list(TopSNP_Ids)[NumberOfSNPS]
    SNPID = ""
    SourceContent = ""
    if "rs" in str(SNP_Ids):
        SNPID = str(SNP_Ids)
        SNPID = SNPID.replace("GSA-", "")
    else:
        try:
            # the below search term is built as per the protocols set by the "Advanced Search Builder tool" in NCBI dbSNP website
            SearchTerm = f"https://www.ncbi.nlm.nih.gov/snp/?term=%28{BP}%5BBase%20Position%20Previous%5D%29%20AND%20{Chr}%5BChromosome%5D"
            # the below line get the HTML source content
            SourceContent = get_source(SearchTerm)
            # the below line takes that source content & tries to filter through to the RS ID of the SNP
            SNPID = get_rsid(SourceContent)
        except Exception:
            # in case the allmighty CGNAT proxy decides to intervene, as usual, as is it's nature, then we retry till we get that source.
            print(f"{str(NumberOfSNPS+1)})Network error : base pair {BP} in chromosome {Chr}\tRetrying.")
            time.sleep(0.1)
            continue
    print(f"{str(NumberOfSNPS+1)})For base pair {BP} in chromosome {Chr}, the RS ID is\t:\t{SNPID}")
    SNPIDList.append(str(SNPID))
    time.sleep(0.1)# better to wait than to anger the Michigan Gods & invite an IP block
    NumberOfSNPS+=1
    if NumberOfSNPS >= HowManyTopSNPIDs:
        break

SNPIDs_RestOfThem = df["SNP"][HowManyTopSNPIDs:]
SNPIDsNew = pd.Series(SNPIDList)
SNPIDs_NewCol = pd.Series(pd.concat([SNPIDsNew,SNPIDs_RestOfThem]))
# SNPIDsNew[0:HowManyTopSNPIDs] = pd.Series(SNPIDList)#[0:len(SNPIDList)+1]

df = df.drop("SNP", 1)
df.insert(1, "SNP", SNPIDs_NewCol)
df.to_csv(f"{inputPath}/{inputFile}_WithRSid.csv", index=False)
