import os

# lhc15o
runlist = [
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431995",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431994",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431993",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431992",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431991",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431989",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431988",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431987",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431986",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431985",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431984",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431983",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431982",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431981",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431980",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431979",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431978",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431977",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431976",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431975",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431974",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431973",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431972",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431971",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431970",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431969",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431968",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431967",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431966",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431965",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431964",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431963",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431962",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431961",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431960",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431957",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431956",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431954",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431953",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431952",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431950",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431949",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431948",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431947",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431946",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431945",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431944",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431943",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431942",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431941",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431940",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431939",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431938",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431937",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431936",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431935",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431934",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431933",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431932",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431931",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431930",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431929",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431928",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431927",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431926",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431925",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431924",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431923",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431922",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431921",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431920",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431919",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431918",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431917",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431916",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431915",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431914",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431913",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431912",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431911",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431910",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431909",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431908",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431907",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431905",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431904",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431903",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431902",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431900",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431899",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431898",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431897",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431896",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431895",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431894",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431893",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431892",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431891",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431890",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431889",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431888",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431887",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431886",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431885",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431884",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431882",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431881",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431880",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431878",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431877",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431876",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431875",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431874",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431873",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431872",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431871",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431870",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431869",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431868",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431867",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431866",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431865",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431864",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431863",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431862",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431861",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431860",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431859",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431858",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431857",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431856",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431855",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431854",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431853",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431852",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431851",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431850",
"/alice/cern.ch/user/a/alihyperloop/jobs/0043/hy_431849",

]

# # lhc18q
# runlist = [
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415029",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415028",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415027",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415025",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415024",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415023",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415022",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415020",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415019",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415015",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415014",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415013",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415012",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415011",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415010",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415009",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415006",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415005",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415004",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415003",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415002",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415001",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415000",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414999",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414998",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414997",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414996",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414995",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414994",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414993",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414992",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414991",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414990",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414989",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414988",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414987",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414986",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414985",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414983",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414982",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414980",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414979",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414978",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414977",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414976",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414975",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414974",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414973",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414971",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414970",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414969",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414968",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414967",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414966",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414964",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414963",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414962",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414961",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414960",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414959",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414957",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414956",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414955",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414954",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414953",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414952",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414951",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414950",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414949",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414948",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414947",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414946",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414945",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414944",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414942",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414941",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414940",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414938",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414937",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414936",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414935",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414934",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414933",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414932",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414931",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414930",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414929",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414928",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414927",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414926",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414925",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414924",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414923",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414922",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414921",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414920",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414919",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414918",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414917",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414915",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414914",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414913",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414912",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414911",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414910",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414909",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414908",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414907",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414906",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414905",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414904",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414903",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414902",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414901",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414900",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414899",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414898",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414897",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414896",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414895",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414894",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414893",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414892",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414891",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_414890"
# ]

# lhc18r
# runlist = [
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415458",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415457",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415456",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415455",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415454",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415453",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415452",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415451",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415450",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415449",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415447",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415446",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415445",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415444",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415443",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415442",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415441",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415440",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415439",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415438",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415437",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415436",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415434",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415433",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415432",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415430",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415429",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415428",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415427",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415426",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415425",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415424",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415423",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415422",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415421",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415420",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415419",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415418",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415417",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415416",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415413",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415412",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415411",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415410",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415409",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415408",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415407",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415406",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415405",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415404",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415403",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415402",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415401",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415400",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415399",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415397",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415396",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415395",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415394",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415393",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415392",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415391",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415390",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415389",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415388",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415387",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415386",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415385",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415384",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415383",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415382",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415381",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415380",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415379",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415378",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415377",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415376",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415375",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415372",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415371",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415370",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415369",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415368",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415366",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415365",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415364",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415362",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415361",
# "/alice/cern.ch/user/a/alihyperloop/jobs/0041/hy_415360"
# ]

for irun, run in enumerate(runlist):
    os.system(f"alien_cp alien://{run}/AnalysisResults.root file:files_15o/AnalysisResults_{irun}.root")