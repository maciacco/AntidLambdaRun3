import os

# lhc15o
# runlist = [
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397532",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397533",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397603",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397543",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397544",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397565",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397576",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397578",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397583",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397594",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397607",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397608",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397611",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397612",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397616",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397537",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397538",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397539",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397540",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397541",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397542",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397548",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397555",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397577",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397584",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397586",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397587",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397588",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397589",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397599",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397602",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397604",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397605",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397628",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397629",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397653",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397654",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397519",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397520",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397530",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397534",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397535",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397536",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397622",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397627",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397630",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397638",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397645",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397652",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397655",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397657",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397658",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397512",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397513",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397529",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397631",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397644",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397514",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397515",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397516",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397521",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397522",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397523",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397524",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397525",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397545",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397546",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397547",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397549",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397550",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397551",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397552",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397553",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397554",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397556",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397557",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397558",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397559",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397560",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397561",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397562",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397563",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397566",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397567",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397579",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397592",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397593",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397595",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397597",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397598",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397600",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397606",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397610",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397613",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397614",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397615",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397620",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397621",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397632",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397633",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397636",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397637",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397639",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397640",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397641",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397648",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397649",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397517",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397569",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397570",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397571",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397572",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397573",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397574",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397575",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397580",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397581",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397582",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397526",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397527",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397528",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397585",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397590",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397601",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397623",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397624",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397625",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397626",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397634",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397635",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397646",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397647",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397650",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397591",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397609",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397617",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397618",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_397619"
# ]

# # lhc18q
runlist = [
    "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398380",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398379",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398378",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398376",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398375",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398374",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398373",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398371",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398370",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398366",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398365",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398364",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398363",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398362",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398361",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398360",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398357",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398356",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398355",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398354",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398353",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398352",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398351",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398350",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398349",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398348",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398347",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398346",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398345",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398344",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398343",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398342",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398341",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398340",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398339",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398338",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398337",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398336",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398334",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398333",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398331",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398330",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398329",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398328",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398327",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398326",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398325",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398324",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398322",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398321",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398320",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398319",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398318",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398317",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398315",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398314",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398313",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398312",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398311",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398310",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398308",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398307",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398306",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398305",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398304",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398303",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398302",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398301",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398300",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398299",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398298",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398297",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398296",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398295",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398293",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398292",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398291",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398289",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398288",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398287",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398286",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398285",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398284",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398283",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398282",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398281",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398280",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398279",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398278",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398277",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398276",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398275",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398274",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398273",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398272",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398271",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398270",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398269",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398268",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398266",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398265",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398264",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398263",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398262",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398261",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398260",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398259",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398258",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398257",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398256",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398255",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398254",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398253",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398252",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398251",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398250",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398249",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398248",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398247",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398246",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398245",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398244",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398243",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398231",
"/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398230",

]

# lhc18r
# runlist = [
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398497",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398496",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398495",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398494",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398493",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398492",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398491",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398490",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398489",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398488",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398486",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398485",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398484",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398483",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398482",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398481",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398480",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398479",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398478",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398477",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398476",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398475",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398473",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398472",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398471",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398469",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398468",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398467",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398466",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398465",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398464",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398463",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398462",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398461",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398460",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398459",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398458",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398457",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398456",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398455",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398452",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398451",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398450",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398449",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398448",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398447",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398446",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398445",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398444",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398443",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398442",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398441",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398440",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398439",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398438",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398436",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398435",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398434",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398433",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398432",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398431",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398430",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398429",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398428",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398427",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398426",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398425",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398424",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398423",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398422",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398421",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398420",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398419",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398418",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398417",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398416",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398415",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398414",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398411",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398410",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398409",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398408",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398407",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398405",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398404",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398403",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398401",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398400",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0039/hy_398399"
# ]

for irun, run in enumerate(runlist):
    os.system(f"alien_cp alien://{run}/AnalysisResults.root file:files_18q/AnalysisResults_{irun}.root")