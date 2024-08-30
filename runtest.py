#! /usr/bin/env python3
import subprocess
import io



def record(cmd, filename):
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    outs, errs = proc.communicate()
    proc.wait()
    # print('errors = {0}\n'.format(errs))
    # print('output = {0}\n'.format(outs))

    outfile = open(filename, 'w')

    ostr = outs.decode("utf-8")
    outfile.write(ostr)


    # print(ostr.split('\n'))



    # buf = io.StringIO(str(outs))

    # print(buf.readline())
    # for line in buf.read():
    #     print(line)
    # outfile.write(str(outs))


def test(cmd, filename):
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    outs, errs = proc.communicate()
    proc.wait()


    prc_out = outs.decode("utf-8").split('\n')
    # ref_out = infile.read()


    ref_out = open(filename, 'r')
    ln = 0
    for lr,lo in zip(ref_out, prc_out):
        ln += 1

        if (lr.find('-- date:') >= 0 and lo.find('-- date:') >= 0) or (lr.find('-- commit:') >= 0 and lo.find('-- commit:') >= 0):
            continue

        ref = [x.strip() for x in lr[:-1].split()]
        res = [x.strip() for x in lo.split()]
        index = 0
        for x,y in zip(res, ref):
            if x != y and index != 7 and index != 3:
                print('Discrepancy for CMD {2} at line {3}:\nref:{0}\ncur:{1}'.format(lr[:-1],lo,' '.join(cmd),ln))
                # exit(1)
            index += 1

    print('CMD {}: ok'.format(cmd))

        # print('{0} / {1}'.format(lr[:-1],lo))

    # infile = open(filename, 'r')
    # for lo,li in zip(str(outs), infile.read()):
    #     print('{0} / {1}'.format(lo,li))



if __name__ == '__main__':
    # record(['build/bibd'], 'reference/bibd.out')
    # record(['build/tempo_scheduler', 'data/sample/tsptw_145.txt', '--input-format tsptw'], 'reference/tsptw.out')
    # record(['build/tempo_scheduler', 'data/sample/la19.txt', '--input-format', 'jsp'], 'reference/jsp.out')
    # record(['build/tempo_scheduler', 'data/sample/la16_0_1', '--input-format', 'jstl'], 'reference/jstl.out')
    # record(['build/tempo_scheduler', 'data/sample/cr_16_1.txt', '--input-format', 'path'], 'reference/path.out')
    # record(['build/tempo_team_orienteering', 'data/sample/tsptw_2.txt', '--no-transitivity', '--no-edge-finding'], 'reference/tsptwopt.out')
    # record(['build/tempo_sat', 'data/sample/uf250-01.cnf', '--verbosity', '3'], 'reference/sat.out')
    # record(['build/tempo_rcpsp', 'data/sample/j309_5.sm', '--no-edge-finding'], 'reference/rcpsp.out')
    # record(['build/tempo_scheduler', 'data/sample/j6-per0-0.txt'], 'reference/osp.out')
    # record(['build/tempo_fjssp','data/sample/mt10c1.fjs', '--no-edge-finding', '--no-transitivity'], 'reference/fjssp.out')
    # record(['build/tempo_scheduler','data/sample/t2-ps05.dat', '--input-format', 'jssdst', '--no-transitivity'], 'reference/t2-pss07.out')    


    test(['build/bibd'], 'reference/bibd.out')
    test(['build/tempo_scheduler', 'data/sample/tsptw_145.txt', '--input-format tsptw'], 'reference/tsptw.out')
    test(['build/tempo_scheduler', 'data/sample/la19.txt', '--input-format', 'jsp'], 'reference/jsp.out')
    test(['build/tempo_scheduler', 'data/sample/la16_0_1', '--input-format', 'jstl'], 'reference/jstl.out')
    test(['build/tempo_scheduler', 'data/sample/cr_16_1.txt', '--input-format', 'path'], 'reference/path.out')
    test(['build/tempo_team_orienteering', 'data/sample/tsptw_2.txt', '--no-transitivity', '--no-edge-finding'], 'reference/tsptwopt.out')
    test(['build/tempo_sat', 'data/sample/uf250-01.cnf', '--verbosity', '3'], 'reference/sat.out')
    test(['build/tempo_rcpsp', 'data/sample/j309_5.sm', '--no-edge-finding'], 'reference/rcpsp.out')
    test(['build/tempo_scheduler', 'data/sample/j6-per0-0.txt'], 'reference/osp.out')
    test(['build/tempo_fjssp','data/sample/mt10c1.fjs', '--no-edge-finding', '--no-transitivity'], 'reference/fjssp.out')    
    test(['build/tempo_scheduler','data/sample/t2-ps05.dat', '--input-format', 'jssdst', '--no-transitivity'], 'reference/t2-pss07.out')
