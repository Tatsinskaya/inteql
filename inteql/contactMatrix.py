class ContactMatrix:
    """
    This class contain all the Hi-C data for a tissue and chromosome. 
    It also, contains the methods to normalize and compute the expected values for contact matrices.
    The information on how to process the data can be found in Rao et. all 2014 supplementary information.
    """
    def __init__(self, matrixFolder, chrm, resl):
        self._chrm = chrm
        self._resl = resl
        self._matrixDirc = matrixFolder
        
        # Observed
        self._matrixPath = self._matrixDirc + '/chr' + str(self._chrm) + '_' + str(resl // 1000) + 'kb.RAWobserved'
        self._krNormPath = self._matrixDirc + '/chr' + str(self._chrm) + '_' + str(resl // 1000) + 'kb.KRnorm'
        self._vcNormPath = self._matrixDirc + '/chr' + str(self._chrm) + '_' + str(resl // 1000) + 'kb.VCnorm'
        self._sqNormPath = self._matrixDirc + '/chr' + str(self._chrm) + '_' + str(resl // 1000) + 'kb.SQRTVCnorm'
        
        # Expected
        self._expectPath = self._matrixDirc + '/chr' + str(self._chrm) + '_' + str(resl // 1000) + 'kb.RAWexpected'
        self._krExpcPath = self._matrixDirc + '/chr' + str(self._chrm) + '_' + str(resl // 1000) + 'kb.KRexpected'
        self._vcExpcPath = self._matrixDirc + '/chr' + str(self._chrm) + '_' + str(resl // 1000) + 'kb.VCexpected'
        self._sqExpcPath = self._matrixDirc + '/chr' + str(self._chrm) + '_' + str(resl // 1000) + 'kb.SQRTVCexpected' 
        

    
    # Hidden functions
    def _readMatrix_toDict(self, matrixPath):
        matrix = {}
        c = 0
        with open(matrixPath, 'r') as f:
            for line in f:
                i, j, m = line.strip().split()
                i, j, m = int(i), int(j), float(m)
                if i not in matrix: matrix[i] = {}
                matrix[i][j] = m
                c+= 1
        print('The raw observed matrix file has %d lines.'%c)
        return matrix
    
    def _readVector_toList(self, vectorPath):
        vector = []
        with open(vectorPath, 'r') as f:
            for line in f:
                vector.append(float(line.strip()))
        return vector
    
    def _normalizeMatrix(self, matrix, vector, r):
        normMatrix = {}
        for i in matrix:
            for j in matrix[i]:
                if i not in normMatrix: normMatrix[i] = {}
                v = i // r
                u = j // r
                normMatrix[i][j] = matrix[i][j] / (vector[v] * vector[u])
        return normMatrix
    
    def _getExpectedMatrix(self, matrix, vector, r):
        expcMatrix = {}
        for i in matrix:
            for j in matrix[i]:
                if i not in expcMatrix: expcMatrix[i] = {}
                l = (j - i) // r
                expcMatrix[i][j] = matrix[i][j] / vector[l]
        return expcMatrix

    def _normalize(self, normType):
        if not hasattr(self, '_rawMatrix'): 
            self._rawMatrix = self._readMatrix_toDict(self._matrixPath)
            
        if normType == 'KR':
            if not hasattr(self, '_krNorm'): self._krNorm = self._readVector_toList(self._krNormPath)
            vector = self._krNorm
                
        elif normType == 'VC':
            if not hasattr(self, '_vcNorm'): self._vcNorm = self._readVector_toList(self._vcNormPath)
            vector = self._vcNorm
        
        elif normType == 'SQ':
            if not hasattr(self, '_sqNorm'): self._sqNorm = self._readVector_toList(self._sqNormPath)
            vector = self._sqNorm
        else: print('Method not known.')
            
        return self._normalizeMatrix(self._rawMatrix, vector, self._resl)
    
    
    def _getExpected(self, expcType):
        if not hasattr(self, '_rawMatrix'):
            self._rawMatrix = self._readMatrix_toDict(self._matrixPath)
            
        if expcType == 'KR':
            if not hasattr(self, '_krExpc'): self._krExpc = self._readVector_toList(self._krExpcPath)
            vector = self._krExpc
                
        elif expcType == 'VC':
            if not hasattr(self, '_vcExpc'): self._vcExpc = self._readVector_toList(self._vcExpcPath)
            vector = self._vcExpc
        
        elif expcType == 'SQ':
            if not hasattr(self, '_sqExpc'): self._sqExpc = self._readVector_toList(self._sqExpcPath)
            vector = self._sqExpc
            
        else: print('Method not known.')
            
        return self._getExpectedMatrix(self._rawMatrix, vector, self._resl)
    
    
    def _normalizeExpected(self, expcType):
        if not hasattr(self, '_rawMatrix'):
            self._rawMatrix = self._readMatrix_toDict(self._matrixPath)
            
        if expcType == 'KR':
            if not hasattr(self, '_krExpc'): self._krExpc = self._readVector_toList(self._krExpcPath)
            if not hasattr(self, '_krNormMatrix'): self._krNormMatrix = self._normalize('KR')
            matrix, vector = self._krNormMatrix, self._krExpc
                
        elif expcType == 'VC':
            if not hasattr(self, '_vcExpc'): self._vcExpc = self._readVector_toList(self._vcExpcPath)
            if not hasattr(self, '_vcNormMatrix'): self._vcNormMatrix = self._normalize('VC')
            matrix, vector = self._vcNormMatrix, self._vcExpc
        
        elif expcType == 'SQ':
            if not hasattr(self, '_sqExpc'): self._sqExpc = self._readVector_toList(self._sqExpcPath)
            if not hasattr(self, '_sqNormMatrix'): self._sqNormMatrix = self._normalize('SQ')
            matrix, vector = self._sqNormMatrix, self._sqExpc
        
        else: print('Method not known.')
            
        return self._getExpectedMatrix(matrix, vector, self._resl)

    
    # Properties
    @property
    def rawMatrix(self):
        if not hasattr(self, '_rawMatrix'):
            self._rawMatrix = self._readMatrix_toDict(self._matrixPath)
        return self._rawMatrix

    @property
    def KRnormMatrix(self):
        if not hasattr(self, '_krNormMatrix'):
            self._krNormMatrix = self._normalize('KR')
        return self._krNormMatrix
    
    @property
    def VCnormMatrix(self):
        if not hasattr(self, '_vcNormMatrix'):
            self._vcNormMatrix = self._normalize('VC')
        return self._vcNormMatrix
    
    @property
    def SQRTVCnormMatrix(self):
        if not hasattr(self, '_sqNormMatrix'):
            self._sqNormMatrix = self._normalize('SQ')
        return self._sqNormMatrix
    
    
    @property
    def KRexpectedMatrix(self):
        if not hasattr(self, '_krExpcMatrix'):
            self._krExpcMatrix = self._getExpected('KR')
        return self._krExpcMatrix
    
    @property
    def VCexpectedMatrix(self):
        if not hasattr(self, '_vcExpcMatrix'):
            self._vcExpcMatrix = self._getExpected('VC')
        return self._vcExpcMatrix
    
    @property
    def SQRTVCexpectedMatrix(self):
        if not hasattr(self, '_sqExpcMatrix'):
            self._sqExpcMatrix = self._getExpected('SQ')
        return self._sqExpcMatrix
    
    @property
    def KRnormExpectedMatrix(self):
        if not hasattr(self, '_krNormExpcMatrix'):
            self._krNormExpcMatrix = self._normalizeExpected('KR')
        return self._krNormExpcMatrix
    
    @property
    def VCnormExpectedMatrix(self):
        if not hasattr(self, '_vcNormExpcMatrix'):
            self._vcNormExpcMatrix = self._normalizeExpected('VC')
        return self._vcNormExpcMatrix
    
    @property
    def SQRTVCnormExpectedMatrix(self):
        if not hasattr(self, '_sqNormExpcMatrix'):
            self._sqNormExpcMatrix = self._normalizeExpected('SQ')
        return self._sqNormExpcMatrix