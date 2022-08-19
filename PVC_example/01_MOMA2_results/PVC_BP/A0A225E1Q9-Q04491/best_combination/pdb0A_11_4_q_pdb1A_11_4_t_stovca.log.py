from pymol.cmd import *
code1 = 'pdb0A_11_4_q'
code2 = 'pdb1A_11_4_t'
cmd.load(code1 + '.pdb' , code1)
cmd.load(code2 + '.pdb' , code2)
cmd.hide('lines', 'all')
cmd.color('blue', code1)
cmd.color('green', code2)
cmd.show('sticks', 'all')
cmd.select('pair_1', '(' + code1 + ' and resi 184) + (' + code2 +' and resi 124)')
cmd.color('red', code1 + ' and resi 184')
cmd.color('orange',  code2 + ' and resi 124')
cmd.select('pair_2', '(' + code1 + ' and resi 185) + (' + code2 +' and resi 125)')
cmd.color('red', code1 + ' and resi 185')
cmd.color('orange',  code2 + ' and resi 125')
cmd.select('pair_3', '(' + code1 + ' and resi 186) + (' + code2 +' and resi 126)')
cmd.color('red', code1 + ' and resi 186')
cmd.color('orange',  code2 + ' and resi 126')
cmd.select('pair_4', '(' + code1 + ' and resi 187) + (' + code2 +' and resi 127)')
cmd.color('red', code1 + ' and resi 187')
cmd.color('orange',  code2 + ' and resi 127')
cmd.select('pair_5', '(' + code1 + ' and resi 188) + (' + code2 +' and resi 128)')
cmd.color('red', code1 + ' and resi 188')
cmd.color('orange',  code2 + ' and resi 128')
cmd.select('pair_6', '(' + code1 + ' and resi 189) + (' + code2 +' and resi 129)')
cmd.color('red', code1 + ' and resi 189')
cmd.color('orange',  code2 + ' and resi 129')
cmd.select('pair_7', '(' + code1 + ' and resi 190) + (' + code2 +' and resi 130)')
cmd.color('red', code1 + ' and resi 190')
cmd.color('orange',  code2 + ' and resi 130')
cmd.select('pair_8', '(' + code1 + ' and resi 193) + (' + code2 +' and resi 136)')
cmd.color('red', code1 + ' and resi 193')
cmd.color('orange',  code2 + ' and resi 136')
cmd.select('pair_9', '(' + code1 + ' and resi 194) + (' + code2 +' and resi 137)')
cmd.color('red', code1 + ' and resi 194')
cmd.color('orange',  code2 + ' and resi 137')
cmd.select('pair_10', '(' + code1 + ' and resi 195) + (' + code2 +' and resi 138)')
cmd.color('red', code1 + ' and resi 195')
cmd.color('orange',  code2 + ' and resi 138')
cmd.select('pair_11', '(' + code1 + ' and resi 196) + (' + code2 +' and resi 139)')
cmd.color('red', code1 + ' and resi 196')
cmd.color('orange',  code2 + ' and resi 139')
cmd.select('pair_12', '(' + code1 + ' and resi 197) + (' + code2 +' and resi 140)')
cmd.color('red', code1 + ' and resi 197')
cmd.color('orange',  code2 + ' and resi 140')
cmd.select('pair_13', '(' + code1 + ' and resi 203) + (' + code2 +' and resi 146)')
cmd.color('red', code1 + ' and resi 203')
cmd.color('orange',  code2 + ' and resi 146')
cmd.select('pair_14', '(' + code1 + ' and resi 204) + (' + code2 +' and resi 147)')
cmd.color('red', code1 + ' and resi 204')
cmd.color('orange',  code2 + ' and resi 147')
cmd.select('pair_15', '(' + code1 + ' and resi 205) + (' + code2 +' and resi 148)')
cmd.color('red', code1 + ' and resi 205')
cmd.color('orange',  code2 + ' and resi 148')
cmd.select('pair_16', '(' + code1 + ' and resi 206) + (' + code2 +' and resi 149)')
cmd.color('red', code1 + ' and resi 206')
cmd.color('orange',  code2 + ' and resi 149')
cmd.select('pair_17', '(' + code1 + ' and resi 207) + (' + code2 +' and resi 150)')
cmd.color('red', code1 + ' and resi 207')
cmd.color('orange',  code2 + ' and resi 150')
cmd.select('pair_18', '(' + code1 + ' and resi 208) + (' + code2 +' and resi 151)')
cmd.color('red', code1 + ' and resi 208')
cmd.color('orange',  code2 + ' and resi 151')
cmd.select('pair_19', '(' + code1 + ' and resi 209) + (' + code2 +' and resi 152)')
cmd.color('red', code1 + ' and resi 209')
cmd.color('orange',  code2 + ' and resi 152')
cmd.select('pair_20', '(' + code1 + ' and resi 210) + (' + code2 +' and resi 153)')
cmd.color('red', code1 + ' and resi 210')
cmd.color('orange',  code2 + ' and resi 153')
cmd.select('pair_21', '(' + code1 + ' and resi 211) + (' + code2 +' and resi 154)')
cmd.color('red', code1 + ' and resi 211')
cmd.color('orange',  code2 + ' and resi 154')
cmd.select('pair_22', '(' + code1 + ' and resi 212) + (' + code2 +' and resi 155)')
cmd.color('red', code1 + ' and resi 212')
cmd.color('orange',  code2 + ' and resi 155')
cmd.select('pair_23', '(' + code1 + ' and resi 213) + (' + code2 +' and resi 156)')
cmd.color('red', code1 + ' and resi 213')
cmd.color('orange',  code2 + ' and resi 156')
cmd.select('pair_24', '(' + code1 + ' and resi 214) + (' + code2 +' and resi 171)')
cmd.color('red', code1 + ' and resi 214')
cmd.color('orange',  code2 + ' and resi 171')
cmd.select('pair_25', '(' + code1 + ' and resi 215) + (' + code2 +' and resi 172)')
cmd.color('red', code1 + ' and resi 215')
cmd.color('orange',  code2 + ' and resi 172')
cmd.select('pair_26', '(' + code1 + ' and resi 216) + (' + code2 +' and resi 173)')
cmd.color('red', code1 + ' and resi 216')
cmd.color('orange',  code2 + ' and resi 173')
cmd.select('pair_27', '(' + code1 + ' and resi 217) + (' + code2 +' and resi 174)')
cmd.color('red', code1 + ' and resi 217')
cmd.color('orange',  code2 + ' and resi 174')
cmd.select('pair_28', '(' + code1 + ' and resi 218) + (' + code2 +' and resi 175)')
cmd.color('red', code1 + ' and resi 218')
cmd.color('orange',  code2 + ' and resi 175')
cmd.select('pair_29', '(' + code1 + ' and resi 219) + (' + code2 +' and resi 176)')
cmd.color('red', code1 + ' and resi 219')
cmd.color('orange',  code2 + ' and resi 176')
cmd.select('pair_30', '(' + code1 + ' and resi 220) + (' + code2 +' and resi 177)')
cmd.color('red', code1 + ' and resi 220')
cmd.color('orange',  code2 + ' and resi 177')
cmd.select('pair_31', '(' + code1 + ' and resi 221) + (' + code2 +' and resi 178)')
cmd.color('red', code1 + ' and resi 221')
cmd.color('orange',  code2 + ' and resi 178')
cmd.select('pair_32', '(' + code1 + ' and resi 222) + (' + code2 +' and resi 179)')
cmd.color('red', code1 + ' and resi 222')
cmd.color('orange',  code2 + ' and resi 179')
cmd.select('pair_33', '(' + code1 + ' and resi 230) + (' + code2 +' and resi 180)')
cmd.color('red', code1 + ' and resi 230')
cmd.color('orange',  code2 + ' and resi 180')
cmd.select('pair_34', '(' + code1 + ' and resi 231) + (' + code2 +' and resi 181)')
cmd.color('red', code1 + ' and resi 231')
cmd.color('orange',  code2 + ' and resi 181')
cmd.select('pair_35', '(' + code1 + ' and resi 232) + (' + code2 +' and resi 182)')
cmd.color('red', code1 + ' and resi 232')
cmd.color('orange',  code2 + ' and resi 182')
cmd.select('pair_36', '(' + code1 + ' and resi 233) + (' + code2 +' and resi 183)')
cmd.color('red', code1 + ' and resi 233')
cmd.color('orange',  code2 + ' and resi 183')
cmd.select('pair_37', '(' + code1 + ' and resi 234) + (' + code2 +' and resi 184)')
cmd.color('red', code1 + ' and resi 234')
cmd.color('orange',  code2 + ' and resi 184')
cmd.select('pair_38', '(' + code1 + ' and resi 235) + (' + code2 +' and resi 185)')
cmd.color('red', code1 + ' and resi 235')
cmd.color('orange',  code2 + ' and resi 185')
cmd.bg_color('gray70')
