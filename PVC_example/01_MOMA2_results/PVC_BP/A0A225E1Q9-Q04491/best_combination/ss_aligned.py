from pymol.cmd import *
cmd.load('pdb0A_1_10_q.pdb', 'pdb0A_1_10_q')
cmd.load('pdb1A_1_10_t.pdb', 'pdb1A_1_10_t')
cmd.color('blue', 'pdb0A_1_10_q')
cmd.color('green', 'pdb1A_1_10_t')
cmd.load('pdb0A_11_4_q.pdb', 'pdb0A_11_4_q')
cmd.load('pdb1A_11_4_t.pdb', 'pdb1A_11_4_t')
cmd.color('blue', 'pdb0A_11_4_q')
cmd.color('green', 'pdb1A_11_4_t')
cmd.load('pdb0A_16_7_q.pdb', 'pdb0A_16_7_q')
cmd.load('pdb1A_16_7_t.pdb', 'pdb1A_16_7_t')
cmd.color('blue', 'pdb0A_16_7_q')
cmd.color('green', 'pdb1A_16_7_t')
cmd.hide('lines', 'all')
cmd.show('cartoon', 'all')
cmd.color('red', 'pdb0A_1_10_q and resi 46-51')
cmd.select('pdb0A_1_10_q_S_1', 'pdb0A_1_10_q and resi 46-51')
cmd.color('red', 'pdb0A_1_10_q and resi 57-61')
cmd.select('pdb0A_1_10_q_S_2', 'pdb0A_1_10_q and resi 57-61')
cmd.color('red', 'pdb0A_1_10_q and resi 75-80')
cmd.select('pdb0A_1_10_q_S_3', 'pdb0A_1_10_q and resi 75-80')
cmd.color('red', 'pdb0A_1_10_q and resi 86-93')
cmd.select('pdb0A_1_10_q_S_4', 'pdb0A_1_10_q and resi 86-93')
cmd.color('red', 'pdb0A_1_10_q and resi 96-102')
cmd.select('pdb0A_1_10_q_S_5', 'pdb0A_1_10_q and resi 96-102')
cmd.color('red', 'pdb0A_1_10_q and resi 105-110')
cmd.select('pdb0A_1_10_q_S_6', 'pdb0A_1_10_q and resi 105-110')
cmd.color('red', 'pdb0A_1_10_q and resi 113-120')
cmd.select('pdb0A_1_10_q_S_7', 'pdb0A_1_10_q and resi 113-120')
cmd.color('red', 'pdb0A_1_10_q and resi 123-131')
cmd.select('pdb0A_1_10_q_S_8', 'pdb0A_1_10_q and resi 123-131')
cmd.color('red', 'pdb0A_1_10_q and resi 147-152')
cmd.select('pdb0A_1_10_q_S_9', 'pdb0A_1_10_q and resi 147-152')
cmd.color('red', 'pdb0A_1_10_q and resi 158-162')
cmd.select('pdb0A_1_10_q_S_10', 'pdb0A_1_10_q and resi 158-162')
cmd.color('orange', 'pdb1A_1_10_t and resi 12-17')
cmd.select('pdb1A_1_10_t_S_1', 'pdb1A_1_10_t and resi 12-17')
cmd.color('orange', 'pdb1A_1_10_t and resi 23-28')
cmd.select('pdb1A_1_10_t_S_2', 'pdb1A_1_10_t and resi 23-28')
cmd.color('orange', 'pdb1A_1_10_t and resi 32-39')
cmd.select('pdb1A_1_10_t_S_3', 'pdb1A_1_10_t and resi 32-39')
cmd.color('orange', 'pdb1A_1_10_t and resi 42-50')
cmd.select('pdb1A_1_10_t_S_4', 'pdb1A_1_10_t and resi 42-50')
cmd.color('orange', 'pdb1A_1_10_t and resi 56-61')
cmd.select('pdb1A_1_10_t_S_5', 'pdb1A_1_10_t and resi 56-61')
cmd.color('orange', 'pdb1A_1_10_t and resi 69-74')
cmd.select('pdb1A_1_10_t_S_6', 'pdb1A_1_10_t and resi 69-74')
cmd.color('orange', 'pdb1A_1_10_t and resi 79-85')
cmd.select('pdb1A_1_10_t_S_7', 'pdb1A_1_10_t and resi 79-85')
cmd.color('orange', 'pdb1A_1_10_t and resi 88-95')
cmd.select('pdb1A_1_10_t_S_8', 'pdb1A_1_10_t and resi 88-95')
cmd.color('orange', 'pdb1A_1_10_t and resi 102-107')
cmd.select('pdb1A_1_10_t_S_9', 'pdb1A_1_10_t and resi 102-107')
cmd.color('orange', 'pdb1A_1_10_t and resi 115-120')
cmd.select('pdb1A_1_10_t_S_10', 'pdb1A_1_10_t and resi 115-120')
cmd.color('red', 'pdb0A_11_4_q and resi 183-188')
cmd.select('pdb0A_11_4_q_S_1', 'pdb0A_11_4_q and resi 183-188')
cmd.color('red', 'pdb0A_11_4_q and resi 196-202')
cmd.select('pdb0A_11_4_q_S_2', 'pdb0A_11_4_q and resi 196-202')
cmd.color('red', 'pdb0A_11_4_q and resi 216-221')
cmd.select('pdb0A_11_4_q_S_3', 'pdb0A_11_4_q and resi 216-221')
cmd.color('red', 'pdb0A_11_4_q and resi 229-235')
cmd.select('pdb0A_11_4_q_S_4', 'pdb0A_11_4_q and resi 229-235')
cmd.color('orange', 'pdb1A_11_4_t and resi 124-129')
cmd.select('pdb1A_11_4_t_S_1', 'pdb1A_11_4_t and resi 124-129')
cmd.color('orange', 'pdb1A_11_4_t and resi 139-142')
cmd.select('pdb1A_11_4_t_S_2', 'pdb1A_11_4_t and resi 139-142')
cmd.color('orange', 'pdb1A_11_4_t and resi 172-177')
cmd.select('pdb1A_11_4_t_S_3', 'pdb1A_11_4_t and resi 172-177')
cmd.color('orange', 'pdb1A_11_4_t and resi 181-188')
cmd.select('pdb1A_11_4_t_S_4', 'pdb1A_11_4_t and resi 181-188')
cmd.color('red', 'pdb0A_16_7_q and resi 266-271')
cmd.select('pdb0A_16_7_q_S_1', 'pdb0A_16_7_q and resi 266-271')
cmd.color('red', 'pdb0A_16_7_q and resi 285-289')
cmd.select('pdb0A_16_7_q_S_2', 'pdb0A_16_7_q and resi 285-289')
cmd.color('red', 'pdb0A_16_7_q and resi 294-303')
cmd.select('pdb0A_16_7_q_S_3', 'pdb0A_16_7_q and resi 294-303')
cmd.color('red', 'pdb0A_16_7_q and resi 306-316')
cmd.select('pdb0A_16_7_q_S_4', 'pdb0A_16_7_q and resi 306-316')
cmd.color('red', 'pdb0A_16_7_q and resi 325-330')
cmd.select('pdb0A_16_7_q_S_5', 'pdb0A_16_7_q and resi 325-330')
cmd.color('red', 'pdb0A_16_7_q and resi 336-341')
cmd.select('pdb0A_16_7_q_S_6', 'pdb0A_16_7_q and resi 336-341')
cmd.color('red', 'pdb0A_16_7_q and resi 355-361')
cmd.select('pdb0A_16_7_q_S_7', 'pdb0A_16_7_q and resi 355-361')
cmd.color('orange', 'pdb1A_16_7_t and resi 207-212')
cmd.select('pdb1A_16_7_t_S_1', 'pdb1A_16_7_t and resi 207-212')
cmd.color('orange', 'pdb1A_16_7_t and resi 220-226')
cmd.select('pdb1A_16_7_t_S_2', 'pdb1A_16_7_t and resi 220-226')
cmd.color('orange', 'pdb1A_16_7_t and resi 231-236')
cmd.select('pdb1A_16_7_t_S_3', 'pdb1A_16_7_t and resi 231-236')
cmd.color('orange', 'pdb1A_16_7_t and resi 244-247')
cmd.select('pdb1A_16_7_t_S_4', 'pdb1A_16_7_t and resi 244-247')
cmd.color('orange', 'pdb1A_16_7_t and resi 257-262')
cmd.select('pdb1A_16_7_t_S_5', 'pdb1A_16_7_t and resi 257-262')
cmd.color('orange', 'pdb1A_16_7_t and resi 269-273')
cmd.select('pdb1A_16_7_t_S_6', 'pdb1A_16_7_t and resi 269-273')
cmd.color('orange', 'pdb1A_16_7_t and resi 278-283')
cmd.select('pdb1A_16_7_t_S_7', 'pdb1A_16_7_t and resi 278-283')
cmd.load('pdb0.pdb', 'pdb0')
cmd.color('blue', 'pdb0')
cmd.hide('lines', 'pdb0')
cmd.remove('solvent')
cmd.show('cartoon', 'pdb0 and chain A')
cmd.center('all')
cmd.zoom('all')
cmd.bg_color('gray70')
cmd.clip('slab',5000)
