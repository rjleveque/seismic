from __future__ import print_function
from pylab import *
from clawpack.pyclaw.gauges import GaugeSolution

outdir1 = '_output'
label1 = 'mapped'
outdir2 = '../water_layer_averaged/_output'
label2 = 'averaged'

figure(310, figsize=(13,8))
clf()
gaugenos = [0,1,2,100,101,102,200,201,202]

for k,gno in enumerate(gaugenos):
    g1 = GaugeSolution(gauge_id=gno, path=outdir1)
    g2 = GaugeSolution(gauge_id=gno, path=outdir2)
    subplot(3,3,k+1)
    plot(g1.t, g1.q[0,:]/1e6, 'r', lw=1.5, label=label1)
    plot(g2.t, g2.q[0,:]/1e6, 'b', lw=0.8, label=label2)
    title('sig11/1e6 at Gauge %s' % gno)
    if k==0: legend(loc='lower left')
    #ylim(-4.5,3)

tight_layout()
savefig('sig11_gauges.png')

figure(311, figsize=(13,8))
clf()
for k,gno in enumerate(gaugenos):
    g1 = GaugeSolution(gauge_id=gno, path=outdir1)
    g2 = GaugeSolution(gauge_id=gno, path=outdir2)
    subplot(3,3,k+1)
    plot(g1.t, g1.q[4,:], 'r', lw=1.5, label=label1)
    plot(g2.t, g2.q[4,:], 'b', lw=0.8, label=label2)
    title('v-velocity at Gauge %s' % gno)
    if k==0: legend(loc='lower left')
    #ylim(-4.5,3)

tight_layout()
savefig('v_gauges.png')
