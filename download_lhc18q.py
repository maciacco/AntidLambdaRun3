import os

# lhc15o
sublist = [
23,
6,
2,
9,
33,
29,
20,
8,
5,
20,
17,
3,
13,
27,
24,
9,
14,
10,
13,
18,
17,
16,
2,
14,
14,
22,
4,
25,
28,
3,
13,
14,
30,
10,
20,
6,
10,
26,
11,
10,
24,
10,
4,
12,
3,
7,
4,
4,
7,
5,
9,
7,
4,
5,
3,
15,
2,
2,
12,
5,
5,
2,
29,
9,
2,
17,
22,
5,
16,
6,
14,
10,
29,
43,
15,
10,
5,
5,
25,
87,
25,
16,
82,
11,
26,
22,
30,
30,
50,
34,
25,
49,
28,
60,
60,
12,
4,
33,
31,
59,
23,
32,
8,
72,
18,
41,
20,
26,
90,
15,
10,
11,
7,
15,
31,
28,
10,
36,
6,
52,
12,
30,
29,
44,
20,
21,
6,
2,
8,
31,
27,
18,
7,
5,
18,
15,
3,
12,
25,
22,
9,
12,
9,
12,
17,
16,
14,
2,
12,
13,
20,
4,
23,
25,
3,
12,
12,
28,
9,
18,
6,
10,
24,
10,
9,
22,
10,
4,
10,
3,
7,
4,
4,
6,
5,
9,
7,
4,
5,
3,
13,
2,
2,
10,
4,
5,
1,
27,
8,
2,
16,
21,
5,
14,
5,
13,
9,
27,
39,
14,
9,
5,
4,
23,
79,
23,
15,
76,
11,
25,
20,
27,
27,
46,
31,
23,
44,
26,
57,
55,
11,
3,
30,
29,
55,
20,
30,
7,
67,
15,
37,
19,
23,
80,
14,
10,
10,
7,
14,
29,
26,
9,
34,
6,
49,
11,
26,
26,
41,
18
]

runlist = [
    "/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435305",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435304",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435303",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435301",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435300",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435299",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435298",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435296",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435295",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435291",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435290",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435289",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435288",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435287",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435286",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435285",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435282",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435281",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435280",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435279",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435278",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435277",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435276",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435275",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435274",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435273",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435272",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435271",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435270",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435269",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435268",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435267",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435266",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435265",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435264",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435263",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435262",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435261",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435259",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435258",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435256",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435255",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435254",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435253",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435252",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435251",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435250",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435249",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435247",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435246",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435245",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435244",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435243",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435242",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435240",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435239",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435238",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435237",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435236",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435235",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435233",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435232",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435231",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435230",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435229",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435228",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435227",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435226",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435225",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435224",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435223",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435222",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435221",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435220",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435218",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435217",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435216",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435214",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435213",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435212",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435211",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435210",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435209",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435208",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435207",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435206",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435205",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435204",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435203",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435202",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435201",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435200",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435199",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435198",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435197",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435196",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435195",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435194",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435193",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435191",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435190",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435189",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435188",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435187",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435186",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435185",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435184",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435183",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435182",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435181",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435180",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435179",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435178",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435177",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435176",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434730",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434729",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434728",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434727",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434726",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434725",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434724",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434723",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434722",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434721",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435437",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435436",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435435",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435433",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435432",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435431",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435430",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435428",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435427",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435423",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435422",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435421",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435420",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435419",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435418",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435417",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435414",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435413",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435412",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435411",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435410",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435409",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435408",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435407",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435406",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435405",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435404",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435403",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435402",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435401",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435400",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435399",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435398",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435397",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435396",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435395",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435394",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435393",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435391",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435390",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435388",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435387",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435386",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435385",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435384",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435383",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435382",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435381",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435379",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435378",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435377",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435376",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435375",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435374",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435372",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435371",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435370",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435369",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435368",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435367",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435365",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435364",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435363",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435362",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435361",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435360",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435359",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435358",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435357",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435356",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435355",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435354",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435353",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435352",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435350",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435349",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435348",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435346",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435345",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435344",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435343",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435342",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435341",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435340",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435339",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435338",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435337",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435336",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435335",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435334",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435333",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435332",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435331",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435330",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435329",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435328",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435327",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435326",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435325",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435323",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435322",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435321",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435320",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435319",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435318",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435317",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435316",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435315",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435314",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435313",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435312",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435311",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435310",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435309",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435308",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435307",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_435306",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434720",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434719",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434718",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434717",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434716",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434715",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434714",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_434713"
]

if (len(runlist) != len(sublist)):
    print("Different sizes")
    exit()

counter = 0
file_cmds = open('cmds_18q', "w")
for irun, run in enumerate(runlist):
    for sub in range(1, sublist[irun] + 1):
        fillz = str(sub).zfill(4)
        counter = counter + 1
        print(f"mkdir -p files_18q/{irun}_{fillz} && alien_cp alien://{run}/{fillz}/AnalysisResults.root file:files_18q/{irun}_{fillz}/AnalysisResults_{irun}_{fillz}.root", file=file_cmds)
print(counter)