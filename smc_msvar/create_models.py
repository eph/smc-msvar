import argparse, numpy as np, pandas as p
from pyvar import MinnesotaPrior, SimsZhaSVARPrior
from fortress import make_smc

#from collections import defaultdict

parser = argparse.ArgumentParser(description=
                                 'Replication of Boganni-Herbst 2017.')

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


if args.prior == 'swz':
    print('You`ve selected the SWZ prior')
    swz_lam = [1,1,1,1.2,0.1,1.0,1.0]
    swz = SimsZhaSVARPrior(dat, swz_lam, p=5, cons=True, presample_moments=premom)

    other_files = {'data.txt': dat,
                   'mu.txt': swz.mu,
                   'sigma.txt': swz.sigma,
                   'qmean.txt': qmean,
                   'qvar.txt': qvar}

    mu_prior = """
       mufile = '{mufile}'
       sigmafile = '{sigmafile}'
       do i = 1, self%ns_mu
          allocate(self%coeff_prior(i)%pr, source=SimsZhaSVARPrior(self%nA+self%nF, mufile, sigmafile))
          npara = npara + self%coeff_prior(i)%pr%npara
       end do
    """.format(mufile='mu.txt', sigmafile='sigma.txt')

    modelfile = open('ms_minnpr.f90','r').read()
    modelfile = modelfile.format(data=dat, prior=args.prior,
                                 m=args.m, v=args.v, p=swz.p,
                                 cons=swz.cons, nA=6, nF=48,
                                 mu_prior=mu_prior)


    lib_path = '/home/eherbst/anaconda3/lib'
    inc_path = '/home/eherbst/anaconda3/include'
    smc = make_smc(modelfile, other_files=other_files, output_directory=args.output_dir,
                   lib_path=lib_path, inc_path=inc_path)


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

    mu_prior = """
          do i = 1,self%ns_mu
             allocate(self%coeff_prior(i)%pr, &
                  source=SVARMinnesotaPrior(self%nobs, self%p, self%constant, &
                  {lam1}, {lam2}, {lam3}, {lam4}, {lam5}, {lamxx}, &
                  {tau}, {ybar}, {sbar}))

             npara = npara + self%coeff_prior(i)%pr%npara
          end do
    """.format(**other_parameters)

    modelfile = open('ms_minnpr.f90','r').read()
    modelfile = modelfile.format(data=dat, prior=args.prior,
                                 m=args.m, v=args.v, p=rfb.p,
                                 cons=rfb.cons, nA=6, nF=rfb.ny**2*rfb.p+rfb.ny,
                                 mufile='mu.txt', sigmafile='sigma.txt',
                                 **other_parameters)
    lib_path = '/home/eherbst/anaconda3/lib'
    include_path = '/home/eherbst/anaconda3/include'

    smc = make_smc(modelfile, other_files=other_files, output_directory=args.output_dir,
                   lib_path=lib_path, inc_path=include_path)
elif args.prior == 'rfb-hier':
    print('You`ve selected the RFB-hier prior') 
    lam = [1.0, 1.2, 1.0, 1.0, 1.0, 1.0]
    rfb = MinnesotaPrior([],lam, presample_moments=premom, p=5)

    other_files = {'data.txt': dat,
                   'qmean.txt': qmean,
                   'qvar.txt': qvar}


    other_parameters = {'lam{:d}'.format(d+1): '{}_wp'.format(val) for d, val in enumerate(lam)}
    other_parameters['tau'] = other_parameters['lam6']
    other_parameters['lamxx'] = '0.1_wp'
    other_parameters['ybar'] = '[-0.0010_wp, 0.0122_wp, 0.0343_wp]' 
    other_parameters['sbar'] = '[ 0.0076_wp, 0.0110_wp, 0.0090_wp]'
    other_parameters['ybar_mean'] = '[0.00_wp, 0.00_wp, 0.00_wp]'
    other_parameters['ybar_std'] = '[0.10_wp, 0.10_wp, 0.10_wp]'

    other_parameters['sbar_s'] = '[0.00658_wp, 0.009526_wp, 0.007794_wp]'
    other_parameters['sbar_nu'] = '[3.00_wp, 3.00_wp, 3.00_wp]'
    other_parameters['hyper_lam_theta'] = '2.61_wp'
    other_parameters['hyper_lam_k'] = '0.61_wp'

    mu_prior = """
          do i = 1,self%ns_mu
             allocate(self%coeff_prior(i)%pr, &
                      source=SVARMinnesotaPriorHyper(self%nobs, self%p, self%constant, &
                      {lam2}, {lam3}, {tau}, {hyper_lam_theta}, {hyper_lam_k},&
                      {hyper_lam_theta}, {hyper_lam_k}, &
                             {ybar_mean}, {ybar_std}, &
             {sbar_s}, {sbar_nu}))

             npara = npara + self%coeff_prior(i)%pr%npara
          end do
    """.format(**other_parameters)

    modelfile = open('ms_minnpr.f90','r').read()
    modelfile = modelfile.format(data=dat, prior=args.prior,
                                 m=args.m, v=args.v, p=rfb.p,
                                 cons=rfb.cons, nA=6, nF=rfb.ny**2*rfb.p+rfb.ny*rfb.cons,
                                 mu_prior=mu_prior)
    lib_path = '/home/eherbst/anaconda3/lib'
    include_path = '/home/eherbst/anaconda3/include'
    smc = make_smc(modelfile, other_files=other_files, output_directory=args.output_dir,
                   lib_path=lib_path, inc_path=include_path)


if args.run:
    results = smc.run(npart=2000, bend=4.0, nphi=2000, nblocks=12,
                      conditional_covariance=True, nproc=4, seed=2601)
    mdd=np.log(results['Z_estimates']).sum()
    print('Estimation of {args.prior:} {args.m:1d}m{args.v:1d}v model complete.\nLog MDD estimate: {mdd:8.3f}'.format(args=args, mdd=mdd))

