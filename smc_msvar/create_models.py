import argparse, numpy as np, pandas as p
from pyvar import MinnesotaPrior, SimsZhaSVARPrior
from fortress import make_smc
from collections import defaultdict

parser = argparse.ArgumentParser(description='Program for estimating MS-VARs via FORTRESS')
parser.add_argument('prior', choices=['swz','rfb','rfb-hier'],
                    help='prior selection')
parser.add_argument('-m', type=int, choices=[1,2], default=1,
                    help='# of mean regimes (default=1)')
parser.add_argument('-v', type=int, choices=[1,2,3,4,5], default=1,
                    help='# of variance regimes (default=1)')
parser.add_argument('--no-run', dest='run', action='store_false',
                    help='don`t run the code (just build the model)')
parser.add_argument('--output-dir', default='_fortress_tmp',
                    help='directory for output')

args = parser.parse_args()


datafile = 'sz_2008_joe_data.csv'
dat = p.read_csv(datafile, names=['xgap','infl','int'])

ybar = [-0.0010, 0.0122, 0.0343]
sbar = [ 0.0076, 0.0110, 0.0090]
premom = np.array([ybar, sbar])

qmean = 5.667*np.eye((args.m))
qmean[qmean==0] = 1

qvar = 5.667*np.eye((args.v))
qvar[qvar==0] = 1

other_parameters = defaultdict(lambda: '-1.0d0')

if args.prior == 'swz':
    print('You`ve selected the SWZ prior')
    swz_lam = [1,1,1,1.2,0.1,1.0,1.0]
    swz = SimsZhaSVARPrior(dat, swz_lam, p=5, cons=True, presample_moments=premom)

    other_files = {'data.txt': dat,
                   'mu.txt': swz.mu,
                   'sigma.txt': swz.sigma,
                   'qmean.txt': qmean,
                   'qvar.txt': qvar}

    [other_parameters[x] for x in
     ['lam1','lam2','lam3','lam4','lam5','lamxx','tau']]
    other_parameters['ybar'] = '[-1.0d0, -1.0d0, -1.0d0]'
    other_parameters['sbar'] = '[-1.0d0, -1.0d0, -1.0d0]'

    modelfile = open('ms_minnpr.f90','r').read()
    modelfile = modelfile.format(data=dat, prior=args.prior,
                                 m=args.m, v=args.v, p=swz.p,
                                 cons=swz.cons, nA=6, nF=48,
                                 mufile='mu.txt', sigmafile='sigma.txt',
                                 **other_parameters)

    
    smc = make_smc(modelfile, other_files=other_files, output_directory=args.output_dir)



elif args.prior == 'rfb':
    print('You`ve selected the RFB prior') 
    lam = [1.0, 1.2, 1.0, 1.0, 1.0, 1.0]
    rfb = MinnesotaPrior([],lam, presample_moments=premom, p=5)
    rfb.ny = 3

    other_files = {'data.txt': dat,
                   'qmean.txt': qmean,
                   'qvar.txt': qvar}


    other_parameters = {'lam{:d}'.format(d+1): '{}_wp'.format(val) for d, val in enumerate(lam)}
    other_parameters['tau'] = other_parameters['lam6']
    other_parameters['lamxx'] = '0.1_wp'
    other_parameters['ybar'] = '[-0.0010_wp, 0.0122_wp, 0.0343_wp]' 
    other_parameters['sbar'] = '[ 0.0076_wp, 0.0110_wp, 0.0090_wp]'

    modelfile = open('ms_minnpr.f90','r').read()
    modelfile = modelfile.format(data=dat, prior=args.prior,
                                 m=args.m, v=args.v, p=rfb.p,
                                 cons=rfb.cons, nA=6, nF=48,
                                 mufile='mu.txt', sigmafile='sigma.txt',
                                 **other_parameters)

    smc = make_smc(modelfile, other_files=other_files, output_directory=args.output_dir)
elif args.prior == 'rfb-hier':
    pass

if args.run:
    results = smc.run(npart=2000, bend=4.0, nphi=2000, nblocks=12,
                      conditional_covariance=True, nproc=4, seed=2601)
    mdd=np.log(results['Z_estimates']).sum()
    print('Estimation of {args.prior:} {args.m:1d}m{args.v:1d}v model complete.\nLog MDD estimate: {mdd:8.3f}'.format(args=args, mdd=mdd))

