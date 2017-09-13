import os

os.system('export FC=gfortran')
libs = ' -lfortress -lfruit -lopenblas -lflap -L/home/eherbst/anaconda3/lib/json-fortran -ljsonfortran'
incs = ' -I/home/eherbst/anaconda3/include/fruit -I/home/eherbst/anaconda3/include/fortress'



names = []
for m,v in [(1,1),(2,1),(1,2)]:
    name = 'swzm{}v{}'.format(m,v)
    od = '{}'.format(name)
    create_str = 'cd .. && python create_models.py swz -m {} -v {} --no-run --output-dir test/{}'
    os.system(create_str.format(m,v,od))

    mod = (open('{}/model_t.f90'.format(od), 'r').read()
           .replace(' model_t', ' model_{}_t'.format(name))
           .replace(' msvar_prior',' msvar_prior{}'.format(name)))

    with open('{}/model_{}_t.f90'.format(od,name), 'w') as f:
        f.write(mod)
    os.system('rm {}/model_t.f90'.format(od))
    os.system('mpif90 -c {}/model_{}_t.f90'.format(od,name)
              + libs + incs + ' -ffree-line-length-1000')
    names.append('model_{}_t.o'.format(name))


from FRUIT import *

test_modules = ['test_swzm1v1.f90']
files = ' test_driver.f90 test_swzm1v1.f90 ' + ' '.join(names) 
suite = test_suite(test_modules)
suite.build_run('test_driver.f90', 'mpif90 -o test_driver '
                + files + libs + incs + ' -ffree-line-length-2000')
suite.summary()
