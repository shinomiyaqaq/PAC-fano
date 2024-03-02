class PCLengthError(Exception):#检查信息位是否超过码长
    """An exception for cases when information_size K
    bigger than codeword length N"""
    pass


class PCLengthDivTwoError(Exception):#检测码长是不是2的次方
    """An exception for cases if codeword length N
     is not divisible by 2"""
    pass


class PCInfoLengthError(Exception):#检查信息位是否设置正确
    """Check if length of information message is equal to
    information_size of the PolarCode object"""