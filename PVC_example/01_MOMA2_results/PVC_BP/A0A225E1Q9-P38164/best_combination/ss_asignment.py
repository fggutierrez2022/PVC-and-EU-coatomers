from pymol import cmd, stored
cmd.load('best_alignment.p1m')
cmd.alter('all', 'ss="L"')
cmd.alter('pdb0A_2_10 and resi 48-52', 'ss = "b"')
cmd.alter('pdb0A_2_10 and resi 56-61', 'ss = "b"')
cmd.alter('pdb0A_2_10 and resi 75-81', 'ss = "b"')
cmd.alter('pdb0A_2_10 and resi 85-92', 'ss = "b"')
cmd.alter('pdb0A_2_10 and resi 98-101', 'ss = "b"')
cmd.alter('pdb0A_2_10 and resi 106-109', 'ss = "b"')
cmd.alter('pdb0A_2_10 and resi 114-119', 'ss = "b"')
cmd.alter('pdb0A_2_10 and resi 149-153', 'ss = "b"')
cmd.alter('pdb0A_2_10 and resi 157-162', 'ss = "b"')
cmd.alter('pdb0A_2_10 and resi 184-189', 'ss = "b"')
cmd.alter('pdb2A_2_10 and resi 56-61', 'ss = "b"')
cmd.alter('pdb2A_2_10 and resi 67-72', 'ss = "b"')
cmd.alter('pdb2A_2_10 and resi 76-81', 'ss = "b"')
cmd.alter('pdb2A_2_10 and resi 141-145', 'ss = "b"')
cmd.alter('pdb2A_2_10 and resi 153-158', 'ss = "b"')
cmd.alter('pdb2A_2_10 and resi 162-167', 'ss = "b"')
cmd.alter('pdb2A_2_10 and resi 177-180', 'ss = "b"')
cmd.alter('pdb2A_2_10 and resi 204-209', 'ss = "b"')
cmd.alter('pdb2A_2_10 and resi 213-218', 'ss = "b"')
cmd.alter('pdb2A_2_10 and resi 222-226', 'ss = "b"')
cmd.alter('pdb0A_12_4 and resi 208-211', 'ss = "b"')
cmd.alter('pdb0A_12_4 and resi 215-222', 'ss = "b"')
cmd.alter('pdb0A_12_4 and resi 229-234', 'ss = "b"')
cmd.alter('pdb0A_12_4 and resi 259-261', 'ss = "b"')
cmd.alter('pdb2A_12_4 and resi 244-247', 'ss = "b"')
cmd.alter('pdb2A_12_4 and resi 253-258', 'ss = "b"')
cmd.alter('pdb2A_12_4 and resi 262-267', 'ss = "b"')
cmd.alter('pdb2A_12_4 and resi 288-290', 'ss = "b"')
cmd.alter('pdb0A_16_7 and resi 266-272', 'ss = "b"')
cmd.alter('pdb0A_16_7 and resi 284-289', 'ss = "b"')
cmd.alter('pdb0A_16_7 and resi 294-302', 'ss = "b"')
cmd.alter('pdb0A_16_7 and resi 307-314', 'ss = "b"')
cmd.alter('pdb0A_16_7 and resi 325-331', 'ss = "b"')
cmd.alter('pdb0A_16_7 and resi 335-342', 'ss = "b"')
cmd.alter('pdb0A_16_7 and resi 354-361', 'ss = "b"')
cmd.alter('pdb2A_16_7 and resi 306-310', 'ss = "b"')
cmd.alter('pdb2A_16_7 and resi 316-321', 'ss = "b"')
cmd.alter('pdb2A_16_7 and resi 325-335', 'ss = "b"')
cmd.alter('pdb2A_16_7 and resi 352-363', 'ss = "b"')
cmd.alter('pdb2A_16_7 and resi 370-377', 'ss = "b"')
cmd.alter('pdb2A_16_7 and resi 381-388', 'ss = "b"')
cmd.alter('pdb2A_16_7 and resi 392-397', 'ss = "b"')
cmd.rebuild()
cmd.set('cartoon_discrete_colors','on')
