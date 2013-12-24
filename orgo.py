import networkx as nx
import xml
import pybel
import openbabel
import PIL
import sys
from PIL import ImageFont
from PIL import Image
from PIL import ImageDraw
from PIL import ImageTk
import time, Tkinter
import copy
class molecule():
	def __init__(self,mol):
		self.react = {}
		self.m = mol
		self.getFGs()
		
	def getFGs(self):
		for fg in fgdict:
			self.react[fg] = []
		self.m.addh()
#		self.m.OBMol.AddHydrogens()
		for fg in fgdict:
			indices = pybel.Smarts(fgdict[fg]).findall(self.m)
			for indset in indices:
				self.react[fg].append([self.m.atoms[ind-1].OBAtom for ind in indset])
		self.m.removeh()
		self.log()
		
	def add_atom(mol,atomicnum):
		newatom = mol.m.OBMol.NewAtom()
		newatom.SetAtomicNum(atomicnum)
		return newatom
	def kill_atom(mol,atom):
		mol.m.OBMol.DeleteAtom(atom)
		mol.getFGs()
	def copy_frag(mol,frag):
		newidx = {}
		for atom in openbabel.OBMolAtomIter(frag.m.OBMol):
			natom = mol.AddAtom(atom.GetAtomicNum())
			newidx[atom.GetIdx()] = natom.GetIdx()
		for bond in openbabel.OBMolBondIter(frag.m.OBMol):
			beg = bond.GetBeginAtomIdx()
			end = bond.GetEndAtomIdx()
			mol.m.OBMol.AddBond(newidx[beg],newidx[end],bond.GetBondOrder())
		return newidx
	
	def add_bond(mol,beginatom,endatom,order=1):
		bond = mol.m.OBMol.GetBond(beginatom.GetIdx(),endatom.GetIdx())
		if bond == None:
			mol.m.OBMol.AddBond(beginatom.GetIdx(),endatom.GetIdx(),order)
		else:
			print "updating"
			bond.SetBondOrder(bond.GetBondOrder()+1)
	def break_bond(mol,beginatom,endatom):
		print beginatom.GetAtomicNum(), endatom.GetAtomicNum()
		bond = mol.m.OBMol.GetBond(beginatom,endatom)
		if bond == None:
			print "Error. No such bond exists"
		else:
			if bond.GetBondOrder() > 1:
				bond.SetBondOrder(bond.GetBondOrder()-1)
			else:
				mol.m.OBMol.DeleteBond(bond)
		
	def connect_frag(mol,othermol,idx1,idx2):
		newidxs = mol.copy_frag(othermol)
		mol.m.OBMol.AddBond(idx1,newidxs[idx2],1)
		return newidxs
		
	def add_from_smile(mol,smile,atom):
		mol.connect_frag(mol_from_smile(smile),atom.GetIdx(),1)
		
	def log(self):
		for atom in openbabel.OBMolAtomIter(self.m.OBMol):
			tstr = "atom " + str(atom.GetIdx()) + " is " + str(atom.GetAtomicNum()) 
			print tstr
		for bond in openbabel.OBMolBondIter(self.m.OBMol):
			tstr = "bond connects " + str(bond.GetBeginAtomIdx()) + " to " + str(bond.GetEndAtomIdx())
			print tstr
		for fg in self.react:
			if len(self.react[fg]) > 0:
				print fg
				i = 0
			for inst in self.react[fg]:
				j = 0
				i = i + 1
				print "Instance" + str(i)
				for ind in inst:
					print "field" + str(j) + "is atom" + str(ind.GetIdx()) + " which is " + str(ind.GetAtomicNum())
					j = j + 1
	def draw(self):	
		self.m.removeh()
		self.m.draw(0,0,1,0)
		# font = ImageFont.truetype("Arial-Bold.ttf",14)
		font = ImageFont.truetype("Arial.ttf",10)
		img=Image.new("RGBA", (500,250),(255,255,255))
		draw = ImageDraw.Draw(img)
		for atom in self.m.atoms:
			if (not atom.type[0] == 'C') or (atom.type[1] == 'l'):
				if atom.type[1].isdigit():
					strg = atom.type[0]
				else:
					strg = atom.type[0] + atom.type[1]
				draw.text((atom.coords[0]*10 + 96.5, atom.coords[1]*10 + 96),strg,(0,0,0),font=font)
			draw.text((atom.coords[0]*10,atom.coords[1]*10),str(atom.OBAtom.GetIdx()),(0,255,0),font=font)
		for bond in openbabel.OBMolBondIter(self.m.OBMol):
			bcoord = self.m.atoms[bond.GetBeginAtomIdx()-1].coords
			ecoord = self.m.atoms[bond.GetEndAtomIdx()-1].coords
			if not self.m.atoms[bond.GetBeginAtomIdx()-1].type[0] == 'C':
				bcoord = [bcoord[0] - .4*(bcoord[0]-ecoord[0]),bcoord[1]-.4*(bcoord[1]-ecoord[1])]
			if not self.m.atoms[bond.GetEndAtomIdx()-1].type[0] == 'C':
				ecoord = [ecoord[0] - .4*(ecoord[0]-bcoord[0]),ecoord[1]-.4*(ecoord[1]-bcoord[1])]
			draw.line([bcoord[0]*10+100,bcoord[1]*10+100,ecoord[0]*10+100,ecoord[1]*10+100],fill='black')
			if (bond.GetBondOrder() > 1):
				newdy = bcoord[0] - ecoord[0]
				newdx = ecoord[1] - bcoord[1]
				draw.line([bcoord[0]*10+100+(newdx*3),bcoord[1]*10+100+(newdy*3),ecoord[0]*10+100+(newdx*3),ecoord[1]*10+100+(newdy*3)],fill='black')
				if bond.GetBondOrder() == 3:
					draw.line([bcoord[0]*10+100-(newdx*3),bcoord[1]*10+100-(newdy*3),ecoord[0]*10+100-(newdx*3),ecoord[1]*10+100-(newdy*3)],fill='black')
		img.show()
		
def mol_from_smile(smile):
	return molecule(pybel.readstring("smi",smile))
	
fgdict = {
"alcohol":"[CX4][OX2;H]",
"la-alcohol":"[OX2][Mg][Br]",
"acyl halide":"[CX3](=[OX1])[F,Cl,Br,I]",
"aldehyde":"[#6][CX3;H1](=O)",
"alkene": "[CX3]=[CX3]",
"alkyne": "[CX2]#[CX2]",
"ahalide": "[CX4][F,Cl,Br,I]",
"amine":"[NX3;H2,H1;!$(NC=O);!$(NS)]",
"cbxacid": "[CX3](=O)[OX2;H1]",
"grignard": "C[Mg][Br]",
"ketone":"[#6][CX3](=O)[#6]",
"nitrile":"C#N",
"phenyl": "c0ccccc0"
}	

def grig_carbonyl(carbonyl,grig,carb_inds,grig_inds):
	grig.break_bond(grig_inds[0],grig_inds[1])
	new_inds = carbonyl.connect_frag(grig,carb_inds[1],grig_inds[0])
	carbonyl.break_bond(carb_inds[1],carb_inds[2])
	carbonyl.add_bond(carb_inds[2],new_inds[grig_inds[1]])
	
def alc_to_carb(alcohol,alcohol_inds):
	if alcohol_inds[0].GetValence() < 4:
		alcohol.add_bond(alcohol_inds[0],alcohol_inds[1])

def quench(grigd,grigd_inds):
	grigd.break_bond(grigd_inds[0],grigd_inds[1])

def brominalcohol(alcohol,alcohol_inds):
	alcohol.break_bond(alcohol_inds[0],alcohol_inds[1])
	idx = alcohol.add_atom(35)
	alcohol.add_bond(alcohol_inds[0],idx)
	alcohol.kill_atom(alcohol_inds[1])
	
def chlorinalcohol(alcohol,alcohol_inds):
	alcohol.break_bond(alcohol_inds[0],alcohol_inds[1])
	idx = alcohol.add_atom(17)
	alcohol.add_bond(alcohol_inds[0],idx)
	alcohol.kill_atom(alcohol_inds[1])

def chlorinacid(acid,acid_inds):
	acid.break_bond(acid_inds[0],acid_inds[2])
	idx = acid.add_atom(17)
	acid.add_bond(acid_inds[0],idx)
	acid.kill_atom(acid_inds[2])
	
def alcoacid(alcohol,alcohol_inds):
	if alcohol_inds[0].GetValence() < 3:
		idx = alcohol.add_atom(8)
		alcohol.add_bond(alcohol_inds[0],idx,order=2)
	else:
		alc_to_carb(alcohol,alcohol_inds)
		
def aldeacid(ald,ald_inds):
	idx = ald.add_atom(8)
	ald.add_bond(ald_inds[1],idx)
	
def acidalc(cbx,cbx_inds):
	cbx.break_bond(cbx_inds[0],cbx_inds[1])
	cbx.break_bond(cbx_inds[0],cbx_inds[1])
	cbx.kill_atom(cbx_inds[1])
	
def carbalc(carb,carb_inds):
	carb.break_bond(carb_inds[1],carb_inds[2])
	
def e1elim(alc,alc_inds):
	max_val = 0
	best = None
	for atom in openbabel.OBAtomAtomIter(alc_inds[0]):
		if atom.GetAtomicNum() == 6 and atom.GetValence() < 4 and atom.GetValence() > max_val:
			max_val = atom.GetValence()
			best = atom
	if max_val == 0:
		print "No beta hydrogens, foo"
	else:
		print max_val, best
		alc.break_bond(alc_inds[0],alc_inds[1])
		alc.add_bond(best,alc_inds[0],2)
		alc.kill_atom(alc_inds[1])
		
def wolfkish(carbonyl,carbonyl_inds):
	carbonyl.break_bond(carbonyl_inds[1],carbonyl_inds[2])
	carbonyl.break_bond(carbonyl_inds[1],carbonyl_inds[2])
	carbonyl.kill_atom(carbonyl_inds[2])
	
def antimarkovoh(alk,alk_inds):
	alk.break_bond(alk_inds[0],alk_inds[1])
	newoh = alk.add_atom(8)
	if (alk_inds[0].GetValence() < alk_inds[1].GetValence()):
		alk.add_bond(alk_inds[0],newoh)
	else:
		alk.add_bond(alk_inds[1],newoh)

def ozonolysis(alk,alk_inds):
	alk.break_bond(alk_inds[0],alk_inds[1])
	alk.break_bond(alk_inds[0],alk_inds[1])
	newoh1 = alk.add_atom(8)
	newoh2 = alk.add_atom(8)
	alk.add_bond(alk_inds[0],newoh1,2)
	alk.add_bond(alk_inds[1],newoh2,2)

def form_grig(ahal,ahal_inds):
	newmg = ahal.add_atom(12)
	ahal.break_bond(ahal_inds[0],ahal_inds[1])
	ahal.add_bond(ahal_inds[0],newmg)
	ahal.add_bond(ahal_inds[1],newmg)
	
reactivity_dict = {"grignard": {"aldehyde": grig_carbonyl, "ketone": grig_carbonyl}}
	
templatedict = {
"Mg0":{"ahalide":form_grig},
"1. O3 2. Me2S":{"alkene":ozonolysis},
"H2NNH2,KOH":{"aldehyde":wolfkish,"ketone":wolfkish},
"1. BH3, THF\n 2. H2O2,NaOH":{"alkene":antimarkovoh},
"DMSO, COCl2, Et3N": {"alcohol": alc_to_carb},
"PBr3": {"alcohol": brominalcohol},
"SOCl2": {"alcohol": chlorinalcohol,"cbxacid": chlorinacid},
"H2CrO4, H2SO4": {"alcohol": alcoacid, "aldehyde": aldeacid},
"conc H2SO4": {"alcohol": e1elim},
"NaBH4": {"aldehyde": carbalc, "ketone": carbalc},
"H+ workup": {"la-alcohol": quench, "grignard": quench},
"LAH": {"aldehyde": carbalc, "ketone": carbalc, "cbxacid": acidalc, "acyl halide": acidalc},
"Wittig": {"carbonyl": wittig}
}

def react_wsmile(sm_mol,added):
	if added in templatedict:
		for possible in templatedict[added]:
			for instance in sm_mol.react[possible]:
				print "HEY"
				templatedict[added][possible](sm_mol,instance)
	else:		
		added_mol = molecule(pybel.readstring("smi",added))
		added_mol.m.removeh()
		for fg in added_mol.react:
			for instance in added_mol.react[fg]:
				for possible in reactivity_dict[fg]:
		#			print sm_mol.react["aldehyde"]
					for instancesm in sm_mol.react[possible]:
						reactivity_dict[fg][possible](sm_mol,added_mol,instancesm,instance)
	sm_mol.getFGs()
def boxproblem_forward(sm,steps):
	for step in steps:
		sm = react_wsmile(sm,step)
	return sm
	
class Application(Tkinter.Frame):
    def mol_load(self,fromline=True):
		if (fromline):
			line = self.e.get()
			self.disp_mol = molecule(pybel.readstring("smi",line))
		self.disp_mol.m.draw(0,"temp.png",1,0)
		time.sleep(.01)
		try:
			image1.close()
		except:
			pass
		image1 = Image.open("temp.png")
		print image1
		tkpi = ImageTk.PhotoImage(image1)
		try:
			label_image.configure(image = tkpi)
		except:
			label_image = Tkinter.Label(self.master, image=tkpi)
		label_image.im = tkpi
		label_image.grid(row=1,column=3,columnspan=6, rowspan=6, sticky=Tkinter.S+Tkinter.E, padx=50, pady=50)
		self.master.title("Image loaded!")
		self.generate_rxns()
		
    def generate_rxns(self):
		for n in self.nb:
			n.destroy()
		self.nb = []
		cur_row = 1
		for mol in templatedict:
			can_react = False
			for possible in templatedict[mol]:
				if (can_react):
					break
				for instance in self.disp_mol.react[possible]:
					newb = Tkinter.Button(self,text=mol,command=lambda mol=mol: self.window_react(self.disp_mol,mol))
					cur_row = cur_row + 1
					self.nb.append(newb)
					self.nb[-1].grid(row=cur_row,column=3)
					can_react = True
					break
		
    def window_react(self,firstmol,smile):
	#	firstmol.log()
		print "reacting with" + smile
		react_wsmile(firstmol,smile)
		self.update()
		
    def update(self):
		self.mol_load(fromline=False)
		
    def createWidgets(self):
        self.QUIT = Tkinter.Button(self)
        self.QUIT["text"] = "QUIT"
        self.QUIT["fg"]   = "red"
        self.QUIT["command"] =  self.quit

        self.QUIT.grid(column=0,row=0)

        self.load = Tkinter.Button(self)
        self.load["text"] = "load mol",
        self.load["command"] = self.mol_load

        self.load.grid(column=1,row=0)
		
        self.e = Tkinter.Entry(self)
        self.e.grid(column=2,row=0)

    def __init__(self, master=None):
        Tkinter.Frame.__init__(self, master)
        self.master = master
        master.geometry("%dx%d+0+0" % (700, 400))
        self.pack()
        self.createWidgets()
        self.nb = []

top = Tkinter.Tk()
root = Tkinter.Toplevel()
top.withdraw()
app = Application(master=root)
app.mainloop()
root.destroy()	
