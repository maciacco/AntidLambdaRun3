import ROOT

ROOT.gRandom.SetSeed(42)

cpa = [0.995, 0.996, 0.997, 0.998]
dcadau = [0.2, 0.3, 0.4, 0.5]
dcav0pv = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
mass = [1, 1.5, 2, 2.5, 3]
dcapv = [0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
ntpccls = [60, 65, 70, 75, 80, 85, 90]
tpcnsig = [1.5, 2, 2.5]
tofnsig = [2.5, 3]
zvtx = [8, 9, 10]

print(f'cpa\tdcadau\tdcav0pv\tmass\tdcapv\tntpccls\ttpcnsig\ttofnsig\tzvtx')
for i in range(30):
    cpa_i = cpa[ROOT.gRandom.Integer(len(cpa))]
    dcadau_i = dcadau[ROOT.gRandom.Integer(len(dcadau))]
    dcav0pv_i = dcav0pv[ROOT.gRandom.Integer(len(dcav0pv))]
    mass_i = mass[ROOT.gRandom.Integer(len(mass))]
    dcapv_i = dcapv[ROOT.gRandom.Integer(len(dcapv))]
    ntpccls_i = ntpccls[ROOT.gRandom.Integer(len(ntpccls))]
    tpcnsig_i = tpcnsig[ROOT.gRandom.Integer(len(tpcnsig))]
    tofnsig_i = tofnsig[ROOT.gRandom.Integer(len(tofnsig))]
    zvtx_i = zvtx[ROOT.gRandom.Integer(len(zvtx))]
    print(f'cpa{i+1} {cpa_i}\tdcadau{i+1} {dcadau_i}\tdcav0pv{i+1} {dcav0pv_i}\tmass{i+1} {mass_i}\tdcapv{i+1} {dcapv_i}\tntpccls{i+1} {ntpccls_i}\ttpcnsig{i+1} {tpcnsig_i}\ttofnsig{i+1} {tofnsig_i}\tzvtx{i+1} {zvtx_i}')
