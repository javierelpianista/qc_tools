def str_to_int(string):
    if isinstance(string, str):
        str_l = string.lower().replace(' ','')

        if 't' in str_l:
            mult = 1024**4
        elif 'g' in str_l:
            mult = 1024**3
        elif 'm' in str_l:
            mult = 1024**2
        elif 'k' in str_l:
            mult = 1024
        else:
            mult = 1

        num = int("".join([s for s in str_l if s.isdigit()]))

    elif isinstance(string, int):
        num = string
        mult = 1

    else: raise Exception('Don\'t know what to do with input')

    return(num*mult)
        
