class atomwiseList(object):
    @staticmethod
    def toGroupOfThree(L):
        """
        Example:

        >>> toGroupOfThree([1,2,3,4,5,6,7,8,9])
        [(1, 2, 3), (4, 5, 6), (7, 8, 9)]

        """

        assert len(L) % 3 == 0
        l = list(zip(L[::3],L[1::3],L[2::3]))

        return l

    @staticmethod
    def toList(L):
        """
        Example:

        >>> toList([(1, 2, 3), (4, 5, 6), (7, 8, 9)])
        [1, 2, 3, 4, 5, 6, 7, 8, 9]

        """
        newL = []
        for group in L:
            for item in group:
                newL.append(item)
        return newL

if __name__ == "__main__":
    import doctest, utils
    doctest.testmod(utils, verbose=False)
