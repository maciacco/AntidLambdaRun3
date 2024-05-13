import os

# lhc15o
sublist = [
29,
23,
7,
33,
26,
23,
24,
13,
14,
57,
35,
40,
20,
52,
33,
19,
8,
18,
12,
27,
26,
51,
13,
9,
31,
23,
32,
8,
16,
3,
29,
32,
13,
3,
33,
9,
20,
32,
9,
4,
54,
6,
8,
71,
17,
39,
55,
8,
79,
18,
41,
5,
18,
17,
23,
77,
14,
14,
21,
84,
21,
54,
75,
28,
35,
7,
8,
16,
76,
34,
31,
20,
20,
8,
15,
7,
11,
44,
7,
10,
12,
18,
31,
67,
21,
9,
79,
14,
29,
23,
18,
6,
26,
20,
18,
19,
11,
12,
47,
28,
33,
16,
42,
26,
16,
6,
14,
9,
21,
21,
40,
10,
7,
24,
19,
26,
8,
12,
3,
23,
25,
11,
3,
26,
7,
16,
25,
8,
3,
43,
6,
6,
57,
13,
30,
45,
7,
62,
15,
33,
4,
15,
13,
18,
62,
12,
11,
16,
68,
17,
43,
60,
23,
28,
6,
5,
14,
60,
27,
25,
16,
16,
8,
12,
6,
9,
35,
6,
8,
10,
15,
25,
53,
17,
7,
62,
12,
24
]

runlist = [
    "/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435528",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435527",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435526",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435525",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435524",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435523",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435522",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435521",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435520",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435519",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435517",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435516",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435515",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435514",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435513",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435512",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435511",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435510",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435509",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435508",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435507",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435506",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435504",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435503",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435502",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435500",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435499",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435498",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435497",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435496",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435495",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435494",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435493",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435492",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435491",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435490",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435489",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435488",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435487",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435486",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435483",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435482",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435481",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435480",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435479",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435478",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435477",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435476",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435475",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435474",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435473",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435472",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435471",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435470",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435469",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435467",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435466",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435465",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435464",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435463",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435462",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435461",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435460",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435459",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435458",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435457",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435456",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435455",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435454",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435453",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435452",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435451",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435450",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435449",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435448",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435447",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435446",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435445",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435442",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435441",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435440",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435439",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435438",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434711",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434710",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434709",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434707",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434706",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434705",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435622",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435621",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435620",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435619",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435618",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435617",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435616",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435615",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435614",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435613",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435611",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435610",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435609",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435608",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435607",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435606",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435605",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435604",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435603",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435602",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435601",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435600",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435598",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435597",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435596",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435594",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435593",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435592",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435591",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435590",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435589",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435588",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435587",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435586",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435585",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435584",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435583",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435582",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435581",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435580",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435577",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435576",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435575",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435574",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435573",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435572",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435571",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435570",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435569",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435568",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435567",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435566",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435565",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435564",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435563",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435561",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435560",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435559",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435558",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435557",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435556",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435555",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435554",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435553",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435552",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435551",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435550",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435549",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435548",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435547",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435546",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435545",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435544",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435543",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435542",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435541",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435540",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435539",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435536",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435535",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435534",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435533",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435532",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435530",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435529",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434704",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434702",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434701",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434700"
]

if (len(runlist) != len(sublist)):
    print("Different sizes")
    exit()

counter = 0
file_cmds = open('cmds_18r', "w")
for irun, run in enumerate(runlist):
    for sub in range(1, sublist[irun] + 1):
        fillz = str(sub).zfill(4)
        counter = counter + 1
        print(f"mkdir -p files_18r/{irun}_{fillz} && alien_cp alien://{run}/{fillz}/AnalysisResults.root file:files_18r/{irun}_{fillz}/AnalysisResults_{irun}_{fillz}.root", file=file_cmds)
print(counter)