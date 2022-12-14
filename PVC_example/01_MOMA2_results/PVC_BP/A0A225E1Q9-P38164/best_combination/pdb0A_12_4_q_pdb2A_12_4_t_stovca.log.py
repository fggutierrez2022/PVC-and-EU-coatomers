from pymol.cmd import *
code1 = 'pdb0A_12_4_q'
code2 = 'pdb2A_12_4_t'
cmd.load(code1 + '.pdb' , code1)
cmd.load(code2 + '.pdb' , code2)
cmd.hide('lines', 'all')
cmd.color('blue', code1)
cmd.color('green', code2)
cmd.show('sticks', 'all')
cmd.select('pair_1', '(' + code1 + ' and resi 216) + (' + code2 +' and resi 252)')
cmd.color('red', code1 + ' and resi 216')
cmd.color('orange',  code2 + ' and resi 252')
cmd.select('pair_2', '(' + code1 + ' and resi 217) + (' + code2 +' and resi 253)')
cmd.color('red', code1 + ' and resi 217')
cmd.color('orange',  code2 + ' and resi 253')
cmd.select('pair_3', '(' + code1 + ' and resi 218) + (' + code2 +' and resi 254)')
cmd.color('red', code1 + ' and resi 218')
cmd.color('orange',  code2 + ' and resi 254')
cmd.select('pair_4', '(' + code1 + ' and resi 219) + (' + code2 +' and resi 255)')
cmd.color('red', code1 + ' and resi 219')
cmd.color('orange',  code2 + ' and resi 255')
cmd.select('pair_5', '(' + code1 + ' and resi 220) + (' + code2 +' and resi 256)')
cmd.color('red', code1 + ' and resi 220')
cmd.color('orange',  code2 + ' and resi 256')
cmd.select('pair_6', '(' + code1 + ' and resi 221) + (' + code2 +' and resi 257)')
cmd.color('red', code1 + ' and resi 221')
cmd.color('orange',  code2 + ' and resi 257')
cmd.select('pair_7', '(' + code1 + ' and resi 222) + (' + code2 +' and resi 258)')
cmd.color('red', code1 + ' and resi 222')
cmd.color('orange',  code2 + ' and resi 258')
cmd.select('pair_8', '(' + code1 + ' and resi 223) + (' + code2 +' and resi 259)')
cmd.color('red', code1 + ' and resi 223')
cmd.color('orange',  code2 + ' and resi 259')
cmd.select('pair_9', '(' + code1 + ' and resi 227) + (' + code2 +' and resi 260)')
cmd.color('red', code1 + ' and resi 227')
cmd.color('orange',  code2 + ' and resi 260')
cmd.select('pair_10', '(' + code1 + ' and resi 228) + (' + code2 +' and resi 261)')
cmd.color('red', code1 + ' and resi 228')
cmd.color('orange',  code2 + ' and resi 261')
cmd.select('pair_11', '(' + code1 + ' and resi 229) + (' + code2 +' and resi 262)')
cmd.color('red', code1 + ' and resi 229')
cmd.color('orange',  code2 + ' and resi 262')
cmd.select('pair_12', '(' + code1 + ' and resi 230) + (' + code2 +' and resi 263)')
cmd.color('red', code1 + ' and resi 230')
cmd.color('orange',  code2 + ' and resi 263')
cmd.select('pair_13', '(' + code1 + ' and resi 231) + (' + code2 +' and resi 264)')
cmd.color('red', code1 + ' and resi 231')
cmd.color('orange',  code2 + ' and resi 264')
cmd.select('pair_14', '(' + code1 + ' and resi 232) + (' + code2 +' and resi 265)')
cmd.color('red', code1 + ' and resi 232')
cmd.color('orange',  code2 + ' and resi 265')
cmd.select('pair_15', '(' + code1 + ' and resi 233) + (' + code2 +' and resi 266)')
cmd.color('red', code1 + ' and resi 233')
cmd.color('orange',  code2 + ' and resi 266')
cmd.select('pair_16', '(' + code1 + ' and resi 234) + (' + code2 +' and resi 267)')
cmd.color('red', code1 + ' and resi 234')
cmd.color('orange',  code2 + ' and resi 267')
cmd.select('pair_17', '(' + code1 + ' and resi 235) + (' + code2 +' and resi 268)')
cmd.color('red', code1 + ' and resi 235')
cmd.color('orange',  code2 + ' and resi 268')
cmd.select('pair_18', '(' + code1 + ' and resi 236) + (' + code2 +' and resi 269)')
cmd.color('red', code1 + ' and resi 236')
cmd.color('orange',  code2 + ' and resi 269')
cmd.select('pair_19', '(' + code1 + ' and resi 253) + (' + code2 +' and resi 281)')
cmd.color('red', code1 + ' and resi 253')
cmd.color('orange',  code2 + ' and resi 281')
cmd.select('pair_20', '(' + code1 + ' and resi 254) + (' + code2 +' and resi 282)')
cmd.color('red', code1 + ' and resi 254')
cmd.color('orange',  code2 + ' and resi 282')
cmd.select('pair_21', '(' + code1 + ' and resi 255) + (' + code2 +' and resi 283)')
cmd.color('red', code1 + ' and resi 255')
cmd.color('orange',  code2 + ' and resi 283')
cmd.select('pair_22', '(' + code1 + ' and resi 256) + (' + code2 +' and resi 284)')
cmd.color('red', code1 + ' and resi 256')
cmd.color('orange',  code2 + ' and resi 284')
cmd.select('pair_23', '(' + code1 + ' and resi 257) + (' + code2 +' and resi 285)')
cmd.color('red', code1 + ' and resi 257')
cmd.color('orange',  code2 + ' and resi 285')
cmd.select('pair_24', '(' + code1 + ' and resi 258) + (' + code2 +' and resi 286)')
cmd.color('red', code1 + ' and resi 258')
cmd.color('orange',  code2 + ' and resi 286')
cmd.select('pair_25', '(' + code1 + ' and resi 259) + (' + code2 +' and resi 287)')
cmd.color('red', code1 + ' and resi 259')
cmd.color('orange',  code2 + ' and resi 287')
cmd.select('pair_26', '(' + code1 + ' and resi 260) + (' + code2 +' and resi 288)')
cmd.color('red', code1 + ' and resi 260')
cmd.color('orange',  code2 + ' and resi 288')
cmd.select('pair_27', '(' + code1 + ' and resi 261) + (' + code2 +' and resi 289)')
cmd.color('red', code1 + ' and resi 261')
cmd.color('orange',  code2 + ' and resi 289')
cmd.bg_color('gray70')
