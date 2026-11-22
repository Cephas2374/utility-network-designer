import sys, wntr
wn = wntr.network.WaterNetworkModel(sys.argv[1])
res = wntr.sim.EpanetSimulator(wn).run_sim(file_prefix=sys.argv[2])
res.node['head'].iloc[0].to_csv(sys.argv[2] + '_head.csv')
res.link['flowrate'].iloc[0].to_csv(sys.argv[2] + '_flow.csv')
