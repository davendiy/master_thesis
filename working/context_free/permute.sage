from trie import TriedDict
from itertools import product
from concurrent.futures import ProcessPoolExecutor
import logging 


x = 'ababxbabxby'
# noinspection PyUnresolvedReferences
F.<a,b,x,y> = FreeGroup()

a_inv, b_inv, x_inv, y_inv = a^(-1), b^(-1), x^(-1), y ^(-1)


relations = {
    (a, x): (x_inv, b), 
    (a, y): (y_inv, b_inv), 
    (a, y_inv): (x, a_inv), 
    (b, x): (y, b_inv), 
    
}


relations_extended = relations.copy()

for el, val in relations.items(): 
    s, t = el 
    t1, s1 = val 
    relations_extended[s1, t^(-1)] = t1^(-1), s

for el, val in relations_extended.copy().items():
    s, t = el 
    t1, s1 = val 
    relations_extended[s1^(-1), t1^(-1)] = t^(-1), s^(-1)

tmp = TriedDict()
for el, v in relations_extended.items(): 
    tmp[el] = v
    
relations = tmp



def gen_word(i1, i2, j1, j2): 
    return [a]*i1 + [b]*i2 + [x]*j1 + [y]*j2 + [a_inv]*i1 + [b_inv]*i2 + [x_inv]*j1 + [x_inv]*j2


def to_tuple(element): 
    res = []
    for el, degree in element.syllables(): 
        res += [el ^ sign(degree)] * abs(degree)
    
    return tuple(res) 
    

def reduce(word, update_relations=False):
    orig_word = word.copy()
    
    while True: 
        succ = False
        for i in range(len(word)):
            pref, left, val = relations.max_prefix(word[i:])
            if len(pref) > 2: 
                print('prefix', len(pref))
            if pref: 
                word = word[:i] + list(val) + left
                succ = True 
                break 

        if not succ: 
            break 
    
    res = x * x^(-1)
    for el in word: 
        res *= el 
    
    if update_relations: 
        relations[tuple(orig_word)] = to_tuple(res)
    return res


def reduce_particular(i1, i2, j1, j2): 
    
    word = gen_word(i1, i2, j1, j2)
    return reduce(word)




logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger('test')

filehandler = logging.FileHandler('context-free.log')
filehandler.setLevel(logging.INFO)
filehandler.setFormatter(logging.Formatter('%(message)s'))
logger.addHandler(filehandler)
logger.info('test')

triv = x * x_inv


def _task(i1, i2, j1, j2): 
    
    word = gen_word(i1, i2, j1, j2)
    res = triv
    for el in word: 
        res *= el 
    if res == triv: 
        return 
        
    if reduce(word) == triv: 
        logger.info(f"{i1}, {i2}, {j1}, {j2}")



for N in range(10, 60):
    
    logger.info(f"Checking for N={N}...")
    with ProcessPoolExecutor(max_workers=32) as p: 
        
        for (i1, i2, j1, j2) in product(list(range(0, N)), repeat=4):
            if N > 10 and (i1 < N-1 and i2 < N-1 and j1 < N-1 and j2 < N-1): 
                continue
            p.submit(_task, i1, i2, j1, j2)

