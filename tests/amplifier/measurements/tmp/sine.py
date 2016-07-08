from __future__ import division
import math

def generate_sine_table(length=1024):    
    raw_table = []  
    for index, item in enumerate((math.sin(2*math.pi*i/length) for i in xrange(length))):
        if math.modf(item)[0] > 0.5:
            value = hex(int(math.ceil((item*0xFF))))
        else:
            value = hex(int(math.floor((item*0xFF))))
        if divmod(index+1, 16)[-1]:
           raw_table.append(hex(int(item*0xFF)))
        else:
           raw_table.append(hex(int(item*0xFF)) + ',\n')
           
    output_table = []    
    for item in (raw_table[j:j+16] for j in xrange(0, len(raw_table), 16)):
        output_table.append(','.join(item))
    
    print ''.join(output_table)

if __name__ == '__main__':
    generate_sine_table(512)
