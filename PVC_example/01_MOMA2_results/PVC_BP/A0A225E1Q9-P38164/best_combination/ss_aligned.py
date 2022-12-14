from pymol.cmd import *
cmd.load('pdb0A_2_10_q.pdb', 'pdb0A_2_10_q')
cmd.load('pdb2A_2_10_t.pdb', 'pdb2A_2_10_t')
cmd.color('blue', 'pdb0A_2_10_q')
cmd.color('green', 'pdb2A_2_10_t')
cmd.load('pdb0A_12_4_q.pdb', 'pdb0A_12_4_q')
cmd.load('pdb2A_12_4_t.pdb', 'pdb2A_12_4_t')
cmd.color('blue', 'pdb0A_12_4_q')
cmd.color('green', 'pdb2A_12_4_t')
cmd.load('pdb0A_16_7_q.pdb', 'pdb0A_16_7_q')
cmd.load('pdb2A_16_7_t.pdb', 'pdb2A_16_7_t')
cmd.color('blue', 'pdb0A_16_7_q')
cmd.color('green', 'pdb2A_16_7_t')
cmd.hide('lines', 'all')
cmd.show('cartoon', 'all')
cmd.color('red', 'pdb0A_2_10_q and resi 48-52')
cmd.select('pdb0A_2_10_q_b_1', 'pdb0A_2_10_q and resi 48-52')
cmd.color('red', 'pdb0A_2_10_q and resi 56-61')
cmd.select('pdb0A_2_10_q_b_2', 'pdb0A_2_10_q and resi 56-61')
cmd.color('red', 'pdb0A_2_10_q and resi 75-81')
cmd.select('pdb0A_2_10_q_b_3', 'pdb0A_2_10_q and resi 75-81')
cmd.color('red', 'pdb0A_2_10_q and resi 85-92')
cmd.select('pdb0A_2_10_q_b_4', 'pdb0A_2_10_q and resi 85-92')
cmd.color('red', 'pdb0A_2_10_q and resi 98-101')
cmd.select('pdb0A_2_10_q_b_5', 'pdb0A_2_10_q and resi 98-101')
cmd.color('red', 'pdb0A_2_10_q and resi 106-109')
cmd.select('pdb0A_2_10_q_b_6', 'pdb0A_2_10_q and resi 106-109')
cmd.color('red', 'pdb0A_2_10_q and resi 114-119')
cmd.select('pdb0A_2_10_q_b_7', 'pdb0A_2_10_q and resi 114-119')
cmd.color('red', 'pdb0A_2_10_q and resi 149-153')
cmd.select('pdb0A_2_10_q_b_8', 'pdb0A_2_10_q and resi 149-153')
cmd.color('red', 'pdb0A_2_10_q and resi 157-162')
cmd.select('pdb0A_2_10_q_b_9', 'pdb0A_2_10_q and resi 157-162')
cmd.color('red', 'pdb0A_2_10_q and resi 184-189')
cmd.select('pdb0A_2_10_q_b_10', 'pdb0A_2_10_q and resi 184-189')
cmd.color('orange', 'pdb2A_2_10_t and resi 56-61')
cmd.select('pdb2A_2_10_t_b_1', 'pdb2A_2_10_t and resi 56-61')
cmd.color('orange', 'pdb2A_2_10_t and resi 67-72')
cmd.select('pdb2A_2_10_t_b_2', 'pdb2A_2_10_t and resi 67-72')
cmd.color('orange', 'pdb2A_2_10_t and resi 76-81')
cmd.select('pdb2A_2_10_t_b_3', 'pdb2A_2_10_t and resi 76-81')
cmd.color('orange', 'pdb2A_2_10_t and resi 141-145')
cmd.select('pdb2A_2_10_t_b_4', 'pdb2A_2_10_t and resi 141-145')
cmd.color('orange', 'pdb2A_2_10_t and resi 153-158')
cmd.select('pdb2A_2_10_t_b_5', 'pdb2A_2_10_t and resi 153-158')
cmd.color('orange', 'pdb2A_2_10_t and resi 162-167')
cmd.select('pdb2A_2_10_t_b_6', 'pdb2A_2_10_t and resi 162-167')
cmd.color('orange', 'pdb2A_2_10_t and resi 177-180')
cmd.select('pdb2A_2_10_t_b_7', 'pdb2A_2_10_t and resi 177-180')
cmd.color('orange', 'pdb2A_2_10_t and resi 204-209')
cmd.select('pdb2A_2_10_t_b_8', 'pdb2A_2_10_t and resi 204-209')
cmd.color('orange', 'pdb2A_2_10_t and resi 213-218')
cmd.select('pdb2A_2_10_t_b_9', 'pdb2A_2_10_t and resi 213-218')
cmd.color('orange', 'pdb2A_2_10_t and resi 222-226')
cmd.select('pdb2A_2_10_t_b_10', 'pdb2A_2_10_t and resi 222-226')
cmd.color('red', 'pdb0A_12_4_q and resi 208-211')
cmd.select('pdb0A_12_4_q_b_1', 'pdb0A_12_4_q and resi 208-211')
cmd.color('red', 'pdb0A_12_4_q and resi 215-222')
cmd.select('pdb0A_12_4_q_b_2', 'pdb0A_12_4_q and resi 215-222')
cmd.color('red', 'pdb0A_12_4_q and resi 229-234')
cmd.select('pdb0A_12_4_q_b_3', 'pdb0A_12_4_q and resi 229-234')
cmd.color('red', 'pdb0A_12_4_q and resi 259-261')
cmd.select('pdb0A_12_4_q_b_4', 'pdb0A_12_4_q and resi 259-261')
cmd.color('orange', 'pdb2A_12_4_t and resi 244-247')
cmd.select('pdb2A_12_4_t_b_1', 'pdb2A_12_4_t and resi 244-247')
cmd.color('orange', 'pdb2A_12_4_t and resi 253-258')
cmd.select('pdb2A_12_4_t_b_2', 'pdb2A_12_4_t and resi 253-258')
cmd.color('orange', 'pdb2A_12_4_t and resi 262-267')
cmd.select('pdb2A_12_4_t_b_3', 'pdb2A_12_4_t and resi 262-267')
cmd.color('orange', 'pdb2A_12_4_t and resi 288-290')
cmd.select('pdb2A_12_4_t_b_4', 'pdb2A_12_4_t and resi 288-290')
cmd.color('red', 'pdb0A_16_7_q and resi 266-272')
cmd.select('pdb0A_16_7_q_b_1', 'pdb0A_16_7_q and resi 266-272')
cmd.color('red', 'pdb0A_16_7_q and resi 284-289')
cmd.select('pdb0A_16_7_q_b_2', 'pdb0A_16_7_q and resi 284-289')
cmd.color('red', 'pdb0A_16_7_q and resi 294-302')
cmd.select('pdb0A_16_7_q_b_3', 'pdb0A_16_7_q and resi 294-302')
cmd.color('red', 'pdb0A_16_7_q and resi 307-314')
cmd.select('pdb0A_16_7_q_b_4', 'pdb0A_16_7_q and resi 307-314')
cmd.color('red', 'pdb0A_16_7_q and resi 325-331')
cmd.select('pdb0A_16_7_q_b_5', 'pdb0A_16_7_q and resi 325-331')
cmd.color('red', 'pdb0A_16_7_q and resi 335-342')
cmd.select('pdb0A_16_7_q_b_6', 'pdb0A_16_7_q and resi 335-342')
cmd.color('red', 'pdb0A_16_7_q and resi 354-361')
cmd.select('pdb0A_16_7_q_b_7', 'pdb0A_16_7_q and resi 354-361')
cmd.color('orange', 'pdb2A_16_7_t and resi 306-310')
cmd.select('pdb2A_16_7_t_b_1', 'pdb2A_16_7_t and resi 306-310')
cmd.color('orange', 'pdb2A_16_7_t and resi 316-321')
cmd.select('pdb2A_16_7_t_b_2', 'pdb2A_16_7_t and resi 316-321')
cmd.color('orange', 'pdb2A_16_7_t and resi 325-335')
cmd.select('pdb2A_16_7_t_b_3', 'pdb2A_16_7_t and resi 325-335')
cmd.color('orange', 'pdb2A_16_7_t and resi 352-363')
cmd.select('pdb2A_16_7_t_b_4', 'pdb2A_16_7_t and resi 352-363')
cmd.color('orange', 'pdb2A_16_7_t and resi 370-377')
cmd.select('pdb2A_16_7_t_b_5', 'pdb2A_16_7_t and resi 370-377')
cmd.color('orange', 'pdb2A_16_7_t and resi 381-388')
cmd.select('pdb2A_16_7_t_b_6', 'pdb2A_16_7_t and resi 381-388')
cmd.color('orange', 'pdb2A_16_7_t and resi 392-397')
cmd.select('pdb2A_16_7_t_b_7', 'pdb2A_16_7_t and resi 392-397')
cmd.load('pdb0.pdb', 'pdb0')
cmd.color('blue', 'pdb0')
cmd.hide('lines', 'pdb0')
cmd.remove('solvent')
cmd.show('cartoon', 'pdb0 and chain A')
cmd.center('all')
cmd.zoom('all')
cmd.bg_color('gray70')
cmd.clip('slab',5000)
