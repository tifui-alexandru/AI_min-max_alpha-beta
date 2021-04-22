import time
import copy
import pygame
import sys


ADANCIME_MAX = 4

culori = {
    'alb': (255, 255, 255), 
    'gri': (192, 192, 192), 
    'verde': (0, 255, 0), 
    'rosu': (255, 0, 0),
    'galben': (255, 255, 0),
    'violet': (153, 0, 153)
}

class Joc:
    """
    Clasa care defineste jocul. Se va schimba de la un joc la altul.
    """

    # constante
    JMIN = None
    JMAX = None

    ZID = '#'
    BOMBA = 'b'
    BOMBA_INACTIVA = 'B'
    PROTECTIE = 'p'
    LIBER = ' '

    directions = [(-1, 0, 'N'), (1, 0, 'S'), (0, -1, 'V'), (0, 1, 'E')]
 
    def __init__(self, k, harta, k_jucatori = {'1': 0, '2': 0}, prot_jucatori = {'1': 0, '2': 0}, bomba_inactiva = {'1': None, '2': None}):
        self.k = k
        self.harta = harta
        self.NR_LINII = len(harta)
        self.NR_COLOANE = len(harta[0])
        self.k_jucatori = k_jucatori
        self.prot_jucatori = prot_jucatori
        self.bomba_inactiva = bomba_inactiva

    def get_bombe(self, x, y):
        # returneaza lista bombelor active care afecteaza pozitia (x, y)
        ans = []

        for linie in range(x, self.NR_LINII):
            if self.harta[linie][y] == Joc.ZID:
                break
            elif self.harta[linie][y] == Joc.BOMBA:
                ans.append((linie, y))

        for linie in range(x - 1, -1, -1):
            if self.harta[linie][y] == Joc.ZID:
                break
            elif self.harta[linie][y] == Joc.BOMBA:
                ans.append((linie, y))

        for coloana in range(y, self.NR_COLOANE):
            if self.harta[x][coloana] == Joc.ZID:
                break
            elif self.harta[x][coloana] == Joc.BOMBA:
                ans.append((x, coloana))

        for coloana in range(y - 1, -1, -1):
            if self.harta[x][coloana] == Joc.ZID:
                break
            elif self.harta[x][coloana] == Joc.BOMBA:
                ans.append((x, coloana))

        return ans

    def get_pos(self, char):
        for linie in range(self.NR_LINII):
            for coloana in range(self.NR_COLOANE):
                if self.harta[linie][coloana] == char:
                    return (linie, coloana)
        return None

    def deseneaza_grid(self, castigator = None, remiza = False):
        pierzator = None
        if castigator:
            pierzator = self.jucator_opus(castigator)
        
        for linie in range(self.NR_LINII):
            culoare = culori['alb']
            for coloana in range(self.NR_COLOANE):
                if self.harta[linie][coloana] == Joc.ZID:
                    culoare = culori['gri']
                elif self.harta[linie][coloana] == castigator:
                    culoare = culori['verde']
                elif self.harta[linie][coloana] == pierzator:
                    culoare = culori['rosu']
                elif self.harta[linie][coloana] == Joc.BOMBA or len(self.get_bombe(linie, coloana) > 0):
                    culoare = culori['galben']

                if remiza and (self.harta[linie][coloana] == castigator or self.harta[linie][coloana] == pierzator):
                    culoare = culori['violet']

                pygame.draw.rect(self.__class__.display, culoare, self.__class__.celuleGrid[linie][coloana]) #alb = (255,255,255)

        if self.harta[linie][coloana] == '1':
            self.__class__.display.blit(self.__class__.img_1,(coloana * (self.__class__.dim_celula+1), linie * (self.__class__.dim_celula+1)))
        elif self.harta[linie][coloana] == '1':
            self.__class__.display.blit(self.__class__.img_2,(coloana * (self.__class__.dim_celula+1), linie * (self.__class__.dim_celula+1)))
        elif self.harta[linie][coloana] == Joc.BOMBA:
            self.__class__.display.blit(self.__class__.img_bomba,(coloana * (self.__class__.dim_celula+1), linie * (self.__class__.dim_celula+1)))
        elif self.harta[linie][coloana] == Joc.BOMBA_INACTIVA:
            if (linie, coloana) == self.bomba_inactiva['1']:
                self.__class__.display.blit(self.__class__.img_bomba_inactiva1,(coloana * (self.__class__.dim_celula+1), linie * (self.__class__.dim_celula+1)))
            else:
                self.__class__.display.blit(self.__class__.img_bomba_inactiva2,(coloana * (self.__class__.dim_celula+1), linie * (self.__class__.dim_celula+1)))
        elif self.harta[linie][coloana] == Joc.PROTECTIE:
             self.__class__.display.blit(self.__class__.img_protectie,(coloana * (self.__class__.dim_celula+1), linie * (self.__class__.dim_celula+1)))

		# pygame.display.flip()
        pygame.display.update()

    @classmethod
    def jucator_opus(cls, jucator):
        return cls.JMAX if jucator == cls.JMIN else cls.JMIN

    @classmethod
    def initializeaza(cls, display, NR_LINII, NR_COLOANE, dim_celula):
        cls.display = display
        cls.dim_celula = dim_celula
        cls.img_1 = pygame.image.load('jucator1.png')
        cls.img_2 = pygame.image.load('jucator2.png')
        cls.img_bomba = pygame.image.load('bomba.png')
        cls.img_bomba_inactiva1 = pygame.image.load('bomba_inactiva1.png')
        cls.img_bomba_inactiva2 = pygame.image.load('bomba_inactiva2.png')
        cls.img_protectie = pygame.image.load('protectie.png')
        cls.celuleGrid = []  # este lista cu patratelele din grid
        for linie in range(NR_LINII):
            for coloana in range(NR_COLOANE):
                patr = pygame.Rect(coloana*(dim_celula+1), linie * (dim_celula+1), dim_celula, dim_celula)
                cls.celuleGrid.append(patr)

    def final(self):
        if self.prot_jucatori['1'] < 0 and self.prot_jucatori['2'] < 0:
            return 'remiza'
        elif self.prot_jucatori['1'] < 0:
            return '2'
        elif self.prot_jucatori['2'] < 0:
            return '1'
        else:
            return None

    def valid_pos(self, x, y):
        if x < 0 or y < 0:
            return False
        if x > self.NR_LINII or y > self.NR_COLOANE:
            return False

        return self.harta[x][y] == Joc.LIBER or self.harta[x][y] == Joc.PROTECTIE

    def valid_mutare(self, jucator, pozitie_noua, pune_bomba, activeaza_bomba):
        # verifica daca mutarea este valida

        # valid pe tabla
        if self.valid_pos(pozitie_noua[0], pozitie_noua[1]) == False:
            return False

        # se poate obtine din pozitia curenta
        (x, y) = self.get_pos(jucator)
        distx = abs(pozitie_noua[0] - x)
        disty = abs(pozitie_noua[1] - y)

        if distx > 1 or disty > 1:
            return False

        # punem fara sa activam
        if pune_bomba == 1 and activeaza_bomba == 0 and self.bomba_inactiva[jucator]:
            return False

        # nu punem, dar ar trebui
        if self.k_jucatori[jucator] + 1 == self.k and pune_bomba == 0:
            return False

    def muta(self, jucator, pozitie_noua, pune_bomba, activeaza_bomba):
        # face mutarea data ca parametru

        if self.valid_mutare(jucator, pozitie_noua, pune_bomba, activeaza_bomba) == False:
            return None

        # am trecut de validari
        new_harta = copy.deepcopy(self.harta)
        new_prot_jucatori = copy.deepcopy(self.prot_jucatori)
        new_k_jucatori = copy.deepcopy(self.k_jucatori)
        new_bomba_inactiva = copy.deepcopy(self.bomba_inactiva)

        # harta
        (x, y) = self.get_pos(jucator)
        new_harta[x][y] = Joc.LIBER
        new_harta[pozitie_noua[0]][pozitie_noua[1]] = jucator

        # protectii
        if self.harta[pozitie_noua[0]][pozitie_noua[1]] == Joc.PROTECTIE:
            new_prot_jucatori[jucator] += 1

        # activez bomba
        if activeaza_bomba == 1:
            (bx, by) = self.bomba_inactiva[jucator]
            new_harta[bx][by] = Joc.BOMBA
            new_bomba_inactiva[jucator] = None

        # punem bomba
        if pune_bomba == 1:
            new_harta[x][y] = Joc.BOMBA_INACTIVA
            new_bomba_inactiva[jucator] = (x, y)
            new_k_jucatori[jucator] = 0
        else:
            new_k_jucatori[jucator] += 1

        ans = Joc(self.k, new_harta, new_k_jucatori, new_prot_jucatori, new_bomba_inactiva)

        # explodam bombele aferente
        bombe = ans.get_bombe(pozitie_noua[0], pozitie_noua[1])
        for b in bombe:
            ans.explodeaza(b[0], b[1])

        return ans

    def explodeaza(self, x, y):
        # explodeaza o bomba recursiv
        if self.harta[x][y] != Joc.BOMBA:
            return

        self.harta[x][y] = Joc.LIBER

        for i in range(len(self.NR_LINII)):
            if self.harta[i][y] == Joc.BOMBA:
                self.explodeaza(i, y)
            elif self.harta[i][y] == '1': 
                self.prot_jucatori['1'] -= 1
            elif self.harta[i][y] == '2': 
                self.prot_jucatori['2'] -= 1 

        for j in range(len(self.NR_COLOANE)):
            if self.harta[x][j] == Joc.BOMBA:
                self.explodeaza(x, j)
            elif self.harta[x][j] == '1':
                self.prot_jucatori['1'] -= 1
            elif self.harta[x][j] == '2':
                self.prot_jucatori['2'] -= 1


    def mutari(self, jucator):
        # intoarce toate mutarile valide
        l_mutari = []
        
        (x, y) = self.get_pos(jucator)
        for d in self.directions:
            for pune in range(2):
                for activeaza in range(2):
                    joc = self.muta(jucator, (x + d[0], y + d[1]), pune, activeaza)
                    if joc:
                        l_mutari.append(joc)
        return l_mutari


    # def estimeaza_scor(self, adancime):
    #     t_final = self.final()
    #     # if (adancime==0):
    #     if t_final == self.__class__.JMAX:
    #         return (self.__class__.scor_maxim+adancime)
    #     elif t_final == self.__class__.JMIN:
    #         return (-self.__class__.scor_maxim-adancime)
    #     elif t_final == 'remiza':
    #         return 0
    #     else:
    #         return (self.linii_deschise(self.__class__.JMAX) - self.linii_deschise(self.__class__.JMIN))

    def sirAfisare(self):
        sir = ""
        for linie in range(self.NR_LINII):
            for coloana in range(self.NR_COLOANE):
                sir += self.harta[linie][coloana]
            sir += "\n"
        return sir

    def __str__(self):
        return self.sirAfisare()

    def __repr__(self):
        return self.sirAfisare()


class Stare:
	"""
	Clasa folosita de algoritmii minimax si alpha-beta
	Are ca proprietate tabla de joc
	Functioneaza cu conditia ca in cadrul clasei Joc sa fie definiti JMIN si JMAX (cei doi jucatori posibili)
	De asemenea cere ca in clasa Joc sa fie definita si o metoda numita mutari() care ofera lista cu configuratiile posibile in urma mutarii unui jucator
	"""

	def __init__(self, tabla_joc, j_curent, adancime, parinte=None, scor=None):
		self.tabla_joc = tabla_joc
		self.j_curent = j_curent

		# adancimea in arborele de stari
		self.adancime = adancime

		# scorul starii (daca e finala) sau al celei mai bune stari-fiice (pentru jucatorul curent)
		self.scor = scor

		# lista de mutari posibile din starea curenta
		self.mutari_posibile = []

		# cea mai buna mutare din lista de mutari posibile pentru jucatorul curent
		self.stare_aleasa = None

	def mutari(self):
		l_mutari = self.tabla_joc.mutari(self.j_curent)
		juc_opus = Joc.jucator_opus(self.j_curent)
		l_stari_mutari = [
			Stare(mutare, juc_opus, self.adancime-1, parinte=self) for mutare in l_mutari]

		return l_stari_mutari

	def __str__(self):
		sir = str(self.tabla_joc) + "(Juc curent:"+self.j_curent+")\n"
		return sir

	def __repr__(self):
		sir = str(self.tabla_joc) + "(Juc curent:"+self.j_curent+")\n"
		return sir


""" Algoritmul MinMax """


def min_max(stare):

	if stare.adancime == 0 or stare.tabla_joc.final():
		stare.scor = stare.tabla_joc.estimeaza_scor(stare.adancime)
		return stare

	# calculez toate mutarile posibile din starea curenta
	stare.mutari_posibile = stare.mutari()

	# aplic algoritmul minimax pe toate mutarile posibile (calculand astfel subarborii lor)
	mutari_scor = [min_max(mutare) for mutare in stare.mutari_posibile]

	if stare.j_curent == Joc.JMAX:
		# daca jucatorul e JMAX aleg starea-fiica cu scorul maxim
		stare.stare_aleasa = max(mutari_scor, key=lambda x: x.scor)
	else:
		# daca jucatorul e JMIN aleg starea-fiica cu scorul minim
		stare.stare_aleasa = min(mutari_scor, key=lambda x: x.scor)
	stare.scor = stare.stare_aleasa.scor
	return stare


def alpha_beta(alpha, beta, stare):
	if stare.adancime == 0 or stare.tabla_joc.final():
		stare.scor = stare.tabla_joc.estimeaza_scor(stare.adancime)
		return stare

	if alpha > beta:
		return stare  # este intr-un interval invalid deci nu o mai procesez

	stare.mutari_posibile = stare.mutari()

	if stare.j_curent == Joc.JMAX:
		scor_curent = float('-inf')

		for mutare in stare.mutari_posibile:
			# calculeaza scorul
			stare_noua = alpha_beta(alpha, beta, mutare)

			if (scor_curent < stare_noua.scor):
				stare.stare_aleasa = stare_noua
				scor_curent = stare_noua.scor
			if(alpha < stare_noua.scor):
				alpha = stare_noua.scor
				if alpha >= beta:
					break

	elif stare.j_curent == Joc.JMIN:
		scor_curent = float('inf')

		for mutare in stare.mutari_posibile:

			stare_noua = alpha_beta(alpha, beta, mutare)

			if (scor_curent > stare_noua.scor):
				stare.stare_aleasa = stare_noua
				scor_curent = stare_noua.scor

			if(beta > stare_noua.scor):
				beta = stare_noua.scor
				if alpha >= beta:
					break
	stare.scor = stare.stare_aleasa.scor

	return stare


def afis_daca_final(stare_curenta):
	final = stare_curenta.tabla_joc.final()
	if(final):
		if (final == "remiza"):
			print("Remiza!")
		else:
			print("A castigat "+final)

		return True

	return False


class Buton:
	def __init__(self, display=None, left=0, top=0, w=0, h=0, culoareFundal=(53, 80, 115), culoareFundalSel=(89, 134, 194), text="", font="arial", fontDimensiune=16, culoareText=(255, 255, 255), valoare=""):
		self.display = display
		self.culoareFundal = culoareFundal
		self.culoareFundalSel = culoareFundalSel
		self.text = text
		self.font = font
		self.w = w
		self.h = h
		self.selectat = False
		self.fontDimensiune = fontDimensiune
		self.culoareText = culoareText
		# creez obiectul font
		fontObj = pygame.font.SysFont(self.font, self.fontDimensiune)
		self.textRandat = fontObj.render(self.text, True, self.culoareText)
		self.dreptunghi = pygame.Rect(left, top, w, h)
		# aici centram textul
		self.dreptunghiText = self.textRandat.get_rect(center=self.dreptunghi.center)
		self.valoare = valoare

	def selecteaza(self, sel):
		self.selectat = sel
		self.deseneaza()

	def selecteazaDupacoord(self, coord):
		if self.dreptunghi.collidepoint(coord):
			self.selecteaza(True)
			return True
		return False

	def updateDreptunghi(self):
		self.dreptunghi.left = self.left
		self.dreptunghi.top = self.top
		self.dreptunghiText = self.textRandat.get_rect(center=self.dreptunghi.center)

	def deseneaza(self):
		culoareF = self.culoareFundalSel if self.selectat else self.culoareFundal
		pygame.draw.rect(self.display, culoareF, self.dreptunghi)
		self.display.blit(self.textRandat, self.dreptunghiText)


class GrupButoane:
	def __init__(self, listaButoane=[], indiceSelectat=0, spatiuButoane=10, left=0, top=0):
		self.listaButoane = listaButoane
		self.indiceSelectat = indiceSelectat
		self.listaButoane[self.indiceSelectat].selectat = True
		self.top = top
		self.left = left
		leftCurent = self.left
		for b in self.listaButoane:
			b.top = self.top
			b.left = leftCurent
			b.updateDreptunghi()
			leftCurent += (spatiuButoane+b.w)

	def selecteazaDupacoord(self, coord):
		for ib, b in enumerate(self.listaButoane):
			if b.selecteazaDupacoord(coord):
				self.listaButoane[self.indiceSelectat].selecteaza(False)
				self.indiceSelectat = ib
				return True
		return False

	def deseneaza(self):
		# atentie, nu face wrap
		for b in self.listaButoane:
			b.deseneaza()

	def getValoare(self):
		return self.listaButoane[self.indiceSelectat].valoare


############# ecran initial ########################
def deseneaza_alegeri(display, tabla_curenta):
	btn_alg = GrupButoane(
		top=30,
		left=30,
		listaButoane=[
			Buton(display=display, w=80, h=30, text="minimax", valoare="minimax"),
			Buton(display=display, w=80, h=30, text="alphabeta", valoare="alphabeta")
                ],
		indiceSelectat=1)
	btn_juc = GrupButoane(
		top=100,
		left=30,
		listaButoane=[
			Buton(display=display, w=35, h=30, text="1", valoare="x"),
			Buton(display=display, w=35, h=30, text="2", valoare="0")
                ],
		indiceSelectat=0)
	ok = Buton(display=display, top=170, left=30, w=40,
	           h=30, text="ok", culoareFundal=(155, 0, 55))
	btn_alg.deseneaza()
	btn_juc.deseneaza()
	ok.deseneaza()
	while True:
		for ev in pygame.event.get():
			if ev.type == pygame.QUIT:
				pygame.quit()
				sys.exit()
			elif ev.type == pygame.MOUSEBUTTONDOWN:
				pos = pygame.mouse.get_pos()
				if not btn_alg.selecteazaDupacoord(pos):
					if not btn_juc.selecteazaDupacoord(pos):
						if ok.selecteazaDupacoord(pos):
							display.fill((0, 0, 0))  # stergere ecran
							tabla_curenta.deseneaza_grid()
							return btn_juc.getValoare(), btn_alg.getValoare()
		pygame.display.update()

def main():
	# setari interf grafica
    pygame.init()
    pygame.display.set_caption("Tifui Ioan Alexandru Omu' cu bombe")

    # determina k si harta initiala
    # k e pe prima linie din fisier, iar harta pe pozitiile urmatoare

    k = None
    harta = None

    with open("harta.txt") as fd:
        temp = fd.read().strip().split('\n')
        k = int(temp[0])
        harta = temp[1:]

    if k is None:
        print("Eroare in citirea fisierului cu harta\n")
        exit(1)
        
    # dimensiunea ferestrei in pixeli
    nl = len(harta)
    nc = len(harta[0])
    w = 50
    ecran = pygame.display.set_mode(size=(nc*(w+1)-1, nl*(w+1)-1))  # N *w+ N-1= N*(w+1)-1
    Joc.initializeaza(ecran, NR_LINII=nl, NR_COLOANE=nc, dim_celula=w)

    # initializare tabla
    tabla_curenta = Joc(k, harta)
    Joc.JMIN, tip_algoritm = deseneaza_alegeri(ecran, tabla_curenta)
    print(Joc.JMIN, tip_algoritm)

    Joc.JMAX = '2' if Joc.JMIN == '1' else '1'

    print("Tabla initiala")
    print(str(tabla_curenta))

    # creare stare initiala
    stare_curenta = Stare(tabla_curenta, '1', ADANCIME_MAX)
    tabla_curenta.deseneaza_grid()
 
    while True:
        while (stare_curenta.j_curent == Joc.JMIN): # e omul la mutare
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    # iesim din program
                    pygame.quit()
                    sys.exit()

                elif event.type == pygame.MOUSEBUTTONDOWN:
                    pos = pygame.mouse.get_pos()  # coordonatele cursorului
                    tabla_curenta = stare_curenta.tabla_joc
                

                activez_bomba = False

                for i in range(nl):
                    for j in range(nc):
                        

                elif event.type == pygame.MOUSEBUTTONDOWN:

                    pos = pygame.mouse.get_pos()  # coordonatele cursorului la momentul clickului

                    for np in range(len(Joc.celuleGrid)):

                        if Joc.celuleGrid[np].collidepoint(pos):
                            # linie=np//Joc.NR_COLOANE
                            coloana = np % Joc.NR_COLOANE
                            ###############################

                            if stare_curenta.tabla_joc.matr[0][coloana] == Joc.LIBER:
                                niv = 0
                                while True:
                                    if niv == Joc.NR_LINII or stare_curenta.tabla_joc.matr[niv][coloana] != Joc.LIBER:
                                        stare_curenta.tabla_joc.matr[niv - 1][coloana] = Joc.JMIN
                                        stare_curenta.tabla_joc.ultima_mutare = (niv-1, coloana)
                                        break
                                    niv += 1

                                # afisarea starii jocului in urma mutarii utilizatorului
                                print("\nTabla dupa mutarea jucatorului")
                                print(str(stare_curenta))

                                stare_curenta.tabla_joc.deseneaza_grid(coloana_marcaj=coloana)
                                # testez daca jocul a ajuns intr-o stare finala
                                # si afisez un mesaj corespunzator in caz ca da
                                if (afis_daca_final(stare_curenta)):
                                    break

                                # S-a realizat o mutare. Schimb jucatorul cu cel opus
                                stare_curenta.j_curent = Joc.jucator_opus(stare_curenta.j_curent)

        # --------------------------------
        else:  # jucatorul e JMAX (calculatorul)
            # Mutare calculator

            # preiau timpul in milisecunde de dinainte de mutare
            t_inainte = int(round(time.time() * 1000))
            if tip_algoritm == 'minimax':
                stare_actualizata = min_max(stare_curenta)
            else:  # tip_algoritm=="alphabeta"
                stare_actualizata = alpha_beta(-500, 500, stare_curenta)
            stare_curenta.tabla_joc = stare_actualizata.stare_aleasa.tabla_joc

            print("Tabla dupa mutarea calculatorului\n"+str(stare_curenta))

            # preiau timpul in milisecunde de dupa mutare
            t_dupa = int(round(time.time() * 1000))
            print("Calculatorul a \"gandit\" timp de " +
                    str(t_dupa-t_inainte)+" milisecunde.")

            stare_curenta.tabla_joc.deseneaza_grid()
            if (afis_daca_final(stare_curenta)):
                break

            # S-a realizat o mutare. Schimb jucatorul cu cel opus
            stare_curenta.j_curent = Joc.jucator_opus(stare_curenta.j_curent)


if __name__ == "__main__":
	main()
	while True:
		for event in pygame.event.get():
			if event.type == pygame.QUIT:
				pygame.quit()
				sys.exit()
