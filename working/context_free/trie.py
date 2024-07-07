#!/usr/bin/env python3
# -*-encoding: utf-8-*-

# created: 15.01.2021
# Excusa. Quod scripsi, scripsi.

# by d.zashkonyi

from collections import deque


# TODO: implement full dict interface
class TriedDict:
    """ Implementation of Dictionary based on trie a.k.a prefix
    tree as a recursive dictionary.

    Parameters
    ----------
    **words : elements for the result dict

    Examples
    --------
    >>> tr_dict = TriedDict(test=2, test2=3)
    >>> tr_dict['hello there'] = 4
    >>> print(tr_dict)
    TriedDict({'t': {'e': {'s': {'t': {'': 2, '2': {'': 3}}}}}, 'h': {'e': {'l': {'l': {'o': {' ': {'t': {'h': {'e': {'r': {'e': {'': 4}}}}}}}}}}}})
    >>> 'hello there' in tr_dict
    True
    >>> 'test3' in tr_dict
    False
    >>> del tr_dict['hello there']
    >>> 'hello there' in tr_dict
    False
    >>> tr_dict.max_prefix('test4')
    ('test', '4', 2)
    >>> list(tr_dict.items())
    [('test', 2), ('test2', 3)]
    >>> list(tr_dict.keys())
    ['test', 'test2']
    >>> list(tr_dict.values())
    [2, 3]
    """
    _end = ''

    def __init__(self, **words):
        root = dict()
        self._root = root
        self._size = 0
        for word, value in words.items():
            self[word] = value

    def __contains__(self, word):
        cur_dict = self._root
        for letter in word:
            if letter not in cur_dict:
                return False
            cur_dict = cur_dict[letter]
        return self._end in cur_dict

    def __setitem__(self, word, value):
        cur_dict = self._root
        for letter in word:
            cur_dict = cur_dict.setdefault(letter, {})
        if self._end not in cur_dict:
            self._size += 1
        cur_dict[self._end] = value

    def __getitem__(self, word):
        cur_dict = self._root
        for letter in word:
            if letter not in cur_dict:
                raise KeyError(f"There is no word '{word}' in the given trie.")
            cur_dict = cur_dict[letter]

        if self._end not in cur_dict:
            raise KeyError(f"There is no word '{word}' in the given trie.")
        return cur_dict[self._end]

    def __delitem__(self, word):
        cur_dict = self._root
        for letter in word:
            if letter not in cur_dict:
                raise KeyError(f"There is no element '{word}' in given trie.")
            cur_dict = cur_dict[letter]
        if self._end not in cur_dict:
            raise KeyError(f"There is no element '{word}' in given trie.")
        else:
            del cur_dict[self._end]
            self._size -= 1

    def max_prefix(self, string: str):
        """ Get prefix of the given string with maximum length that
        contains in the TriedDict.

        Parameters
        ----------
        string  : str of any length

        Returns
        -------
        (<found_prefix>, <left_string>, <value>)
        """
        cur_dict = self._root
        last_i = 0
        last_value = None
        i = 0
        for i, letter in enumerate(string):
            if self._end in cur_dict:
                last_i = i
                last_value = cur_dict[self._end]
            if letter in cur_dict:
                cur_dict = cur_dict[letter]
            else:
                break
        else:
            if self._end in cur_dict:
                last_i = i + 1
                last_value = cur_dict[self._end]

        return string[:last_i], string[last_i:], last_value

    def __repr__(self):
        return f'TriedDict({self._root})'

    def __len__(self):
        return self._size

    def items(self):
        stack = deque()
        res_word = deque([''])

        stack.append((iter(self._root), self._root))
        while stack:
            cur_root_iter, cur_root = stack.pop()
            try:
                letter = next(cur_root_iter)
                if letter == self._end:
                    yield ''.join(res_word), cur_root[letter]
                    stack.append((cur_root_iter, cur_root))
                else:
                    stack.append((cur_root_iter, cur_root))
                    stack.append((iter(cur_root[letter]), cur_root[letter]))
                    res_word.append(letter)
            except StopIteration:
                res_word.pop()

    def keys(self):
        for key, _ in self.items():
            yield key

    def values(self):
        for _, value in self.items():
            yield value

    def clear(self):
        self._root.clear()

    def __iter__(self):
        return self.keys()


class _deleted_attr:

    def __set_name__(self, owner, name):
        self._name = name

    def __get__(self, instance, owner):
        raise AttributeError(f"'{owner}' object has no attribute '{self._name}'")


class Trie(TriedDict):
    """ Implementation of Trie a.k.a prefix tree as a recursive dictionary.

    Parameters
    ----------
    *words : elements for the resuld trie

    Examples
    --------
    >>> trie = Trie("word", "word2", "word3", "w", "wo", "next")
    >>> print(trie)
    Trie({'w': {'o': {'r': {'d': {'': '', '2': {'': ''}, '3': {'': ''}}}, '': ''}, '': ''}, 'n': {'e': {'x': {'t': {'': ''}}}}})
    >>> 'word' in trie
    True
    >>> 'word123' in trie
    False
    >>> 'wor' in trie
    False
    >>> trie.add('word1231')
    >>> 'word1231' in trie
    True
    >>> trie.remove('w')
    >>> 'w' in trie
    False
    >>> list(trie)
    ['word', 'word2', 'word3', 'word1231', 'wo', 'next']
    >>> trie.max_prefix('word123next')
    ('word', '123next')

    Raises
    ------
    TypeError  : when called methods __getitem__, __setitem__, __delitem__

    Notes
    -----
    Trie is based on TriedDict, however it doesn't support container
    protocol doesn't have attributes keys, values and items.
    """

    def __init__(self, *words):
        super(Trie, self).__init__()
        for word in words:
            self.add(word)

    def add(self, word):
        super(Trie, self).__setitem__(word, self._end)

    def remove(self, word):
        super(Trie, self).__delitem__(word)

    def __setitem__(self, key, value):
        raise TypeError("'Trie' object does not support item assignment" )

    def __getitem__(self, item):
        raise TypeError("'Trie' object is not subscriptable")

    def __delitem__(self, key):
        raise TypeError("'Trie' object does not support item deletion")

    def max_prefix(self, string):
        res_prefix, left_string, _ = super(Trie, self).max_prefix(string)
        return res_prefix, left_string

    def __repr__(self):
        return f'Trie({self._root})'

    def __iter__(self):
        for key, _ in super(Trie, self).items():
            yield key

    keys = _deleted_attr()
    values = _deleted_attr()
    items = _deleted_attr()
