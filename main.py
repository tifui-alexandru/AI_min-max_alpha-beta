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

tip_estimare = 1

human_move = False

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
    BOMBA_EXPLODATA = 'E'
    PROTECTIE = 'p'
    LIBER = ' '

    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)] # N S V E
 
    def __init__(self, k, harta, k_jucatori = {'1': 0, '2': 0}, prot_jucatori = {'1': 0, '2': 0}, bomba_inactiva = {'1': None, '2': None}):
        self.k = k
        self.harta = harta
        self.NR_LINII = len(harta)
        self.NR_COLOANE = len(harta[0])
        self.k_jucatori = k_jucatori
        self.prot_jucatori = prot_jucatori
        self.bomba_inactiva = bomba_inactiva

    def get_bombe(self, x, y):
        # returneaza lista bombelor active care ar fi activate de o mutare la pozitia (x, y)
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

    def patrat_bombat(self, x, y):
        # intoarce True daca celula (x, y) este afectata de o bomba
        # sau daca pe celula (x, y) se afla o bomba
        for i in range(self.NR_LINII):
            if self.harta[i][y] == Joc.BOMBA:
                return True
        for j in range(self.NR_COLOANE):
            if self.harta[x][j] == Joc.BOMBA:
                return True
        return False

    def deseneaza_grid(self, castigator = None, remiza = False):
        pierzator = None
        if castigator:
            pierzator = self.jucator_opus(castigator)
        
        for linie in range(self.NR_LINII):
            for coloana in range(self.NR_COLOANE):
                
                culoare = culori['alb']

                if self.harta[linie][coloana] == Joc.ZID:
                    culoare = culori['gri']
                elif self.harta[linie][coloana] == castigator:
                    culoare = culori['verde']
                elif self.harta[linie][coloana] == pierzator:
                    culoare = culori['rosu']
                elif self.harta[linie][coloana] == Joc.BOMBA_EXPLODATA:
                    culoare = culori['rosu']
                elif self.patrat_bombat(linie, coloana):
                    culoare = culori['galben']

                if remiza and (self.harta[linie][coloana] == '1' or self.harta[linie][coloana] == '2'):
                    culoare = culori['violet']

                pygame.draw.rect(self.__class__.display, culoare, self.__class__.celuleGrid[linie][coloana]) #alb = (255,255,255)

                if self.harta[linie][coloana] == '1':
                    self.__class__.display.blit(self.__class__.img_1,(coloana * (self.__class__.dim_celula+1), linie * (self.__class__.dim_celula+1)))
                elif self.harta[linie][coloana] == '2':
                    self.__class__.display.blit(self.__class__.img_2,(coloana * (self.__class__.dim_celula+1), linie * (self.__class__.dim_celula+1)))
                elif self.harta[linie][coloana] == Joc.BOMBA:
                    self.__class__.display.blit(self.__class__.img_bomba,(coloana * (self.__class__.dim_celula+1), linie * (self.__class__.dim_celula+1)))
                elif self.harta[linie][coloana] == Joc.BOMBA_EXPLODATA:
                    self.__class__.display.blit(self.__class__.img_bomba,(coloana * (self.__class__.dim_celula+1), linie * (self.__class__.dim_celula+1)))
                    if self.final() is None:
                        self.harta[linie][coloana] = Joc.LIBER
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
        cls.img_1 = pygame.transform.scale(cls.img_1, (dim_celula, dim_celula))
        cls.img_2 = pygame.image.load('jucator2.png')
        cls.img_2 = pygame.transform.scale(cls.img_2, (dim_celula, dim_celula))
        cls.img_bomba = pygame.image.load('bomba.png')
        cls.img_bomba = pygame.transform.scale(cls.img_bomba, (dim_celula, dim_celula))
        cls.img_bomba_inactiva1 = pygame.image.load('bomba_inactiva1.png')
        cls.img_bomba_inactiva1 = pygame.transform.scale(
            cls.img_bomba_inactiva1, (dim_celula, dim_celula))
        cls.img_bomba_inactiva2 = pygame.image.load('bomba_inactiva2.png')
        cls.img_bomba_inactiva2 = pygame.transform.scale(
            cls.img_bomba_inactiva2, (dim_celula, dim_celula))
        cls.img_protectie = pygame.image.load('protectie.png')
        cls.img_protectie = pygame.transform.scale(
            cls.img_protectie, (dim_celula, dim_celula))
        cls.celuleGrid = [[None for _ in range(NR_COLOANE)] for _ in range(NR_LINII)]  # este lista cu patratelele din grid
        for linie in range(NR_LINII):
            for coloana in range(NR_COLOANE):
                patr = pygame.Rect(coloana*(dim_celula+1), linie * (dim_celula+1), dim_celula, dim_celula)
                cls.celuleGrid[linie][coloana] = patr

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
        if x >= self.NR_LINII or y >= self.NR_COLOANE:
            return False

        return self.harta[x][y] == Joc.LIBER or self.harta[x][y] == Joc.PROTECTIE

    def valid_mutare(self, jucator, pozitie_noua, pune_bomba, activeaza_bomba):
        # verifica daca mutarea este valida

        global human_move

        human_move = True

        # valid pe tabla
        if self.valid_pos(pozitie_noua[0], pozitie_noua[1]) == False:
            if human_move:
                print("\nPozitie proasta")
            return False

        # se poate obtine din pozitia curenta
        (x, y) = self.get_pos(jucator)
        distx = abs(pozitie_noua[0] - x)
        disty = abs(pozitie_noua[1] - y)

        if distx + disty > 1:
            if human_move:
                print("\nPozitie nevecina")
            return False

        # vreau sa activez, dar nu este
        if activeaza_bomba and self.bomba_inactiva[jucator] is None:
            return False

        # punem fara sa activam
        if pune_bomba == 1 and activeaza_bomba == 0 and self.bomba_inactiva[jucator]:
            if human_move:
                print("\nTrebuia sa activam bomba")
            return False

        # nu punem, dar ar trebui
        if self.k_jucatori[jucator] + 1 == self.k and pune_bomba == 0:
            if human_move:
                print("\nTrebuia sa punem bomba")
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

        self.harta[x][y] = Joc.BOMBA_EXPLODATA

        for i in range(self.NR_LINII):
            if self.harta[i][y] == Joc.BOMBA:
                self.explodeaza(i, y)
            elif self.harta[i][y] == '1': 
                self.prot_jucatori['1'] -= 1
            elif self.harta[i][y] == '2': 
                self.prot_jucatori['2'] -= 1 

        for j in range(self.NR_COLOANE):
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

    def estimare_scor1(self, adancime, jucator):
        '''
        Definim o functie f(patrat), unde patrat este un patrat vecin jucatorului curent
        Aceasta ne va da un scor de 'cat de bine' ar fi sa mutam in patratul respectiv
        Estimarea va fi data de suma acestor scoruri

        Facem notatiile:
             X - un patrat care ar declansa cel putin o bomba 
             Y - un patrat din afara hartii
             NP - numarul de protectii pe care le detine jucatorul curent
             NB - numarul de bombe care s-ar declansa si ar afecta jucatorul de pe patratul X
             proc = 0 daca la pozitia vecina nu se afla o pozitie si 1 altfel

        Atribuim functia in felul urmator:
            f(LIBER) = 100
            f(ZID/BOMBA) = -100
            f(PROTECTIE) = 130
            f(Y) = -100
                    | -100 daca NP + proc < NB
            f(X) =  |  
                    |  100 * (NP + proc - NB) / NP daca NP + proc >= NB

        Estimarea ordoneaza starile deoarece fiecare jucator tinde sa aiba o zona cat mai mare de deplasare.
        De asemena, jucatorii cu mai multe protectii sunt avantajati, deoarece isi extind oarecum zona de deplasare.
        ''' 

        final = self.final()
        if final:
            if final == jucator:
                return float('inf')
            elif final == self.jucator_opus(jucator):
                return float('-inf')
            else:
                return 0

        # altfel adancimea este 0
        (x, y) = self.get_pos(jucator)
        total = 0

        for d in self.directions:
            new_x = x + d[0]
            new_y = y + d[1]
            if new_x < 0 or y < new_y or new_x >= self.NR_LINII or new_y >= self.NR_COLOANE:
                # Y
                total -= 100
            elif self.harta[new_x][new_y] == Joc.BOMBA or self.harta[new_x][new_y] == Joc.ZID or self.harta[new_x][new_y] == Joc.BOMBA_INACTIVA:
                total -= 100
            elif self.patrat_bombat(new_x, new_y):
                proc = 0
                if self.harta[new_x][new_y] == Joc.PROTECTIE:
                    proc = 1

                NP = proc + self.prot_jucatori[jucator]
                NB = 0

                for i in range(self.NR_LINII):
                    if self.harta[i][new_y] == Joc.BOMBA:
                        NB += 1

                for j in range(self.NR_COLOANE):
                    if self.harta[new_x][j] == Joc.BOMBA:
                        NB += 1

                if NP < NB:
                    total -= 100
                else:
                    total += 100 * (NP - NB) / NP
                
            else:
                if self.harta[new_x][new_y] == Joc.LIBER:
                    total += 100
                else:
                    # protectie
                    total += 130

        return total

    def corespunde(self, x, y, tip_campuri):
        # functie pentru pentru
        if x < 0 or y < 0 or x >= self.NR_LINII or y >= self.NR_COLOANE:
            return False
        if self.harta[x][y] not in tip_campuri:
            return False
        if self.patrat_bombat(x, y):
            return False
        return True

    def calculeaza_aria(self, x, y, k, tip_campuri, campuri_de_numarat):
        # functie ajutatoare pentru estimare2
        # vom face un bfs
        vizitat = {(x, y)}
        q = [(x, y, k)] 

        total = 0

        while len(q):
            (x, y, k) = q.pop(0)
            if self.corespunde(x, y, tip_campuri):
                if self.harta[x][y] in campuri_de_numarat:
                    total += 1

                if k == 0:
                    continue
            
                for d in self.directions:
                    new_x = x + d[0]
                    new_y = y + d[0]
                    if (new_x, new_y) not in vizitat:
                        vizitat.add((new_x, new_y))
                        q.append((new_x, new_y, k - 1))
        return total

    def estimare_scor2(self, adancime, jucator):
        '''
        Fie k' = min(k - k_jucator1, k - k_jucator2)
        A = numarul de patrate accesibile de jucator aflate la o distanta <= k'
        NP = numarul de protectii ale jucatorului in momentul actual
        NPA = numarul de protectii aflate la o distanta <= k' de jucator

        Estimarea va fi (A / (k' * (k' + 1)) - 1) * 100 + NP + NPA / 2

        Cum am ajuns la aceasta formula?
        
        Fiecare jucator doreste sa aiba o arie cat mai mare in care sa se poata deplasa
        pana cand o bomba va fi activata, adica in k' mutari.

        De asemenea, acesta este avantajat daca are un numar de protectii.
        Am considerat ca 2 protectii aflate in aria posibila de miscari valoreaza
        cat o protectie deja aflata in posesia jucatorului (practic am aproximat ca probabilitatea
        ca acesta sa ia o protectie este de 1/2).

        Dorim o estimare care pentru A = 0 sa dea -100, iar pentru A = 2 * k'(k' + 1) sa dea 100.
        Observam ca 2 * k'(k' + 1) este valoarea maxima posibila a lui A.
        De aici vine partea de (A' / (k' * (k' + 1)) - 1) * 100 a formulei,
        protectiile reprezentand un bonus.
        '''

        final = self.final()
        if final:
            if final == jucator:
                return float('inf')
            elif final == self.jucator_opus(jucator):
                return float('-inf')
            else:
                return 0

        # altfel adancimea este 0
        k_prim = min(self.k - self.k_jucatori['1'], self.k - self.k_jucatori['2'])
        (x, y) = self.get_pos(jucator)
        A = self.calculeaza_aria(x, y, k_prim, {Joc.LIBER, Joc.PROTECTIE}, {Joc.LIBER, Joc.PROTECTIE})
        NP = self.prot_jucatori[jucator]
        NPA = self.calculeaza_aria(x, y, k_prim, {Joc.LIBER, Joc.PROTECTIE}, {Joc.PROTECTIE})
        
        return (A / (k_prim * (k_prim + 1)) - 1) * 100 + NP + NPA / 2

    def estimeaza_scor(self, adancime):
        # ne intereseaza ca JMAX sa aiba scor mare, iar JMIN scor mic
        # vom returna diferenta celor doua
        global tip_estimare
        if tip_estimare == 1:
            a = self.estimare_scor1(adancime, Joc.JMAX)
            b = self.estimare_scor1(adancime, Joc.JMIN)
        else:
            a = self.estimare_scor2(adancime, Joc.JMAX)
            b = self.estimare_scor2(adancime, Joc.JMIN)
        return a - b

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

        # scorul starii curente
        self.scor = scor

        # lista de mutari posibile din starea curenta
        self.mutari_posibile = []

        # cea mai buna mutare din lista de mutari posibile pentru jucatorul curent
        self.stare_aleasa = None

    def mutari(self):
        l_mutari = self.tabla_joc.mutari(self.j_curent)
        print(len(l_mutari))
        juc_opus = Joc.jucator_opus(self.j_curent)
        l_stari_mutari = [Stare(mutare, juc_opus, self.adancime-1, parinte=self) for mutare in l_mutari]

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


""" Algoritmul AlphaBeta """

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
            stare_curenta.tabla_joc.deseneaza_grid(remiza = True)
            print("Remiza!")
        else:
            stare_curenta.tabla_joc.deseneaza_grid(castigator = final)
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
        indiceSelectat=0)
    btn_juc = GrupButoane(
        top=100,
        left=30,
        listaButoane=[
            Buton(display=display, w=80, h=30, text="Jucator1", valoare="1"),
            Buton(display=display, w=80, h=30, text="Jucator2", valoare="2")
                ],
        indiceSelectat=0)
    btn_mod_joc = GrupButoane(
        top=170,
        left=30,
        listaButoane=[
            Buton(display=display, w=100, h=30, text="OM vs PC", valoare="OM vs PC"),
            Buton(display=display, w=100, h=30, text="OM vs OM", valoare="OM vs OM"),
            Buton(display=display, w=100, h=30, text="PC vs PC", valoare="PC vs PC")
        ],
        indiceSelectat=0)
    btn_estimare = GrupButoane(
        top=240,
        left=30,
        listaButoane=[
            Buton(display=display, w=130, h=30, text="Estimare 1", valoare=1),
            Buton(display=display, w=130, h=30, text="Estimare 2", valoare=2)
                ],
        indiceSelectat=0)
    btn_adancime = GrupButoane(
        top=310,
        left=30,
        listaButoane=[
            Buton(display=display, w=130, h=30, text="Nivel incepator", valoare="incepator"),
            Buton(display=display, w=130, h=30, text="Nivel mediu", valoare="mediu"),
            Buton(display=display, w=130, h=30, text="Nivel avansat", valoare="avansat")
        ],
        indiceSelectat=0)
    
    ok = Buton(display=display, top=400, left=30, w=40,
                h=30, text="ok", culoareFundal=(155, 0, 55))
    btn_alg.deseneaza()
    btn_juc.deseneaza()
    btn_mod_joc.deseneaza()
    btn_estimare.deseneaza()
    btn_adancime.deseneaza()
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
                        if not btn_mod_joc.selecteazaDupacoord(pos):
                            if not btn_estimare.selecteazaDupacoord(pos):
                                if not btn_adancime.selecteazaDupacoord(pos):
                                    if ok.selecteazaDupacoord(pos):
                                        display.fill((0, 0, 0))  # stergere ecran
                                        tabla_curenta.deseneaza_grid()
                                        return (btn_juc.getValoare(), btn_alg.getValoare(), btn_mod_joc.getValoare(), btn_estimare.getValoare(), btn_adancime.getValoare())
        pygame.display.update()

def eval_conditie(mod_joc, stare_curenta):
    # functie care determina cine este la mutare
    if mod_joc == "OM vs OM":
        return True
    if mod_joc == "PC vs PC":
        return False
    return stare_curenta.j_curent == Joc.JMIN

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
        harta = [list(line) for line in temp[1:]]

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
    global tip_estimare
    global ADANCIME_MAX
    (Joc.JMIN, tip_algoritm, mod_joc, tip_estimare, dificultate) = deseneaza_alegeri(ecran, tabla_curenta)

    print("\n\n")
    
    print("Mod de joc:", mod_joc)

    if mod_joc == "OM vs PC":
        print("Jucatorul uman este jucatorul", Joc.JMIN)

    if mod_joc == "OM vs PC" or "PC vs PC":
        print("Algoritm folosit:", tip_algoritm)
        print("Dificultate: nivel ", dificultate)

        if dificultate == "incepator":
            ADANCIME_MAX = 2
        elif dificultate == "mediu":
            ADANCIME_MAX = 3
        else:
            ADANCIME_MAX = 4

        print("Adancime maxima: ", ADANCIME_MAX)

    if mod_joc == "OM vs PC":
        print("Tip estimare: ", tip_estimare)
    elif mod_joc == "PC vs PC":
        print("PC1: estimare1")
        print("PC2: estimare2")

    print("\n")

    Joc.JMAX = '2' if Joc.JMIN == '1' else '1'

    print("Tabla initiala")
    print(str(tabla_curenta))

    # creare stare initiala
    stare_curenta = Stare(tabla_curenta, '1', ADANCIME_MAX)
    tabla_curenta.deseneaza_grid()

    activez_bomba = False
    pun_bomba = False
    pozitie_noua = (0, 0)
    terminat_mutare = False

    while True:
        if eval_conditie(mod_joc, stare_curenta): # e omul la mutare
            global human_move
            human_move = True
            terminat_mutare = False
            while terminat_mutare == False:
                for event in pygame.event.get():
                    if event.type == pygame.QUIT:
                        # iesim din program
                        pygame.quit()
                        sys.exit()

                    elif event.type == pygame.MOUSEBUTTONDOWN:
                        pos = pygame.mouse.get_pos()  # coordonatele cursorului
                        tabla_curenta = stare_curenta.tabla_joc

                        de_activat_bomba = False
                        for i in range(nl):
                            for j in range(nc):
                                if Joc.celuleGrid[i][j].collidepoint(pos):
                                    if tabla_curenta.harta[i][j] == Joc.BOMBA_INACTIVA:
                                        de_activat_bomba = True
                                    pozitie_noua = (i, j)
                        
                        if event.button == 3 or stare_curenta.tabla_joc.k_jucatori[stare_curenta.j_curent] + 1 == k:
                            pun_bomba = True

                        if de_activat_bomba:
                            activez_bomba = True
                            tabla_intermediara = copy.deepcopy(tabla_curenta)
                            tabla_intermediara.harta[pozitie_noua[0]][pozitie_noua[1]] = Joc.BOMBA
                            tabla_intermediara.deseneaza_grid()
                            print("Tabla intermediara")
                            print(str(tabla_intermediara))
                        else:
                            if pun_bomba and stare_curenta.tabla_joc.bomba_inactiva[stare_curenta.j_curent]:
                                activez_bomba = True

                            tabla_noua = stare_curenta.tabla_joc.muta(stare_curenta.j_curent, pozitie_noua, pun_bomba, activez_bomba)
                            
                            activez_bomba = False

                            if tabla_noua is None:
                                print("Mutare invalida\n")
                                continue

                            stare_curenta = Stare(tabla_noua, Joc.jucator_opus(stare_curenta.j_curent), ADANCIME_MAX)

                            tabla_noua.deseneaza_grid()
                            print("Muta jucatorul", Joc.jucator_opus(stare_curenta.j_curent))
                            print("Tabla dupa mutarea omului")
                            print(str(tabla_noua))

                            terminat_mutare = True

                # resetam valorile
                if terminat_mutare:
                    if afis_daca_final(stare_curenta):
                        while True:
                            for event in pygame.event.get():
                                if event.type == pygame.QUIT:
                                    # iesim din program
                                    pygame.quit()
                                    sys.exit()

                    activez_bomba = False
                    pun_bomba = False
                    pozitie_noua = (0, 0)
                    human_move = False
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

            print("Muta jucatorul", Joc.jucator_opus(stare_curenta.j_curent))
            print("Tabla dupa mutarea calculatorului\n"+str(stare_curenta))

            # preiau timpul in milisecunde de dupa mutare
            t_dupa = int(round(time.time() * 1000))
            print("Calculatorul a \"gandit\" timp de " +
                    str(t_dupa-t_inainte)+" milisecunde.")

            stare_curenta.tabla_joc.deseneaza_grid()

            if afis_daca_final(stare_curenta):
                while True:
                    for event in pygame.event.get():
                        if event.type == pygame.QUIT:
                            # iesim din program
                            pygame.quit()
                            sys.exit()

            # S-a realizat o mutare. Schimb jucatorul cu cel opus
            stare_curenta.j_curent = Joc.jucator_opus(stare_curenta.j_curent)

            # schimb tipul de estimare pt urmatoarea mutare daca joaca PC vs PC
            if mod_joc == "PC vs PC":
                tip_estimare = 3 - tip_estimare


if __name__ == "__main__":
	main()
	while True:
		for event in pygame.event.get():
			if event.type == pygame.QUIT:
				pygame.quit()
				sys.exit()
