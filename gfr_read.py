import pydicom as pyd
import pprint


def read_dcm():
#    paths = '/run/user/1000/gvfs/smb-share:server=rghheapac001.regionh.top.local,share=data/TEMP/GFR_test/DICOM/000010B0/AA5E46BF/AA3453A6/00007800/EE7C9DEF' 
#    paths = '/run/user/1000/gvfs/smb-share:server=rghheapac001.regionh.top.local,share=data/TEMP/GFR_test/DICOM/000010B0/AA55D4DE/AAD50535/0000EAEA/EEF30F12' 
    paths = '/run/user/1000/gvfs/smb-share:server=rghheapac001.regionh.top.local,share=data/TEMP/GFR_test/DICOM/000010B0/AA6051AB/AA0BEF99/0000CE00/EECD8F0A' 
#    paths = '/run/user/1000/gvfs/smb-share:server=rghheapac001.regionh.top.local,share=data/TEMP/GFR_test/DICOM/000010B0/AA453894/AADED158/0000449E/EE9C1A05' 



    dataset = pyd.dcmread(paths)
    return dataset

# Keys 0x00231002 (GFR version) and 0x00231040 (Clearance Comment) failed
# Keys 0x00231021 (Sample Time), x00231022 (Count Per Minuts) and 0x00231023 (Deviation on Sample)
# lies in the sequence of Keys 0x0023103F (Clearance History) and 0x00231020 (Clearance Tests)
# Key 0x0023103F (Clearance History) exist in GFR fler punkt, having 0x00231012,0x00231014,
# 0x00231018 and 0x00231020, with 0x00231020 having x00231021,0x00231022,0x00231023

# information that haven't been found in any GFR data yet:
# ('GFR Version',0x00231002), ('Clearance Comment',0x00231040)

# names and keys used in the outermost layer
GFR_keys = [('GFR',0x00231001), ('GFR Version',0x00231002),('GFR Method',0x00231010),('Body Surface Method',0x00231011),
    ('Clearance',0x00231012),('Normalized Clearance',0x00231014),('Injection time',0x00231018),
    ('Vial number',0x00231019),('Injection weight',0x0023101A),('Vial weight before injection',0x0023101B),
    ('Vial weight after injection',0x0023101C),('Clearance Tests',0x00231020),
    ('Standard Counts Per',0x00231024),('Thining Factor',0x00231028),
    ('Examination Status',0x00231032),('Clearance History',0x0023103F), ('Clearance Comment',0x00231040)]

# names and keys used by 0x00231020
key_1020 = [('Sample Time',0x00231021),('Count Per Minuts',0x00231022),
            ('Deviation on Sample',0x00231023)]

# names and keys used by 0x0023103f
key_103f = [('Clearance',0x00231012),('Normalized Clearance',0x00231014),
            ('Injection time',0x00231018),('Clearance Tests',0x00231020)]



def get_gfrdata(dataset):
    """
    -----------
    get_gfrdata
    -----------
    Inputs:
        dataset: DICOM file

    A function which goes through the gfr keys and accesses them.
    In case of extra layer keys (1020 and 103f) it will also go through their layers.
    Currently just prints.
    """

    tuples = []
    for (name, key) in GFR_keys:
        if key in dataset:
            data = dataset[key].value

            if key == 0x00231020:
                data = data[0]
                tupkey1 = []
                for (name1, key1) in key_1020:
                    if key1 in data:
                        tupkey1.append((name1, data[key1].value))
                    else:
                        print(name1 +' is not in ' + name)
                tuples.append((name, dict(tupkey1)))

            elif key == 0x0023103f:
                data = data[0]
                tupkey1 = []
                for (name1, key1) in key_103f:
                    if key1 in data:
                        data1 = data[key1].value
                        if key1 == 0x00231020:
                            tupkey2 = []
                            data1 = data1[0]
                            for (name2, key2) in key_1020:
                                if key2 in data1:
                                    tupkey2.append((name2, data1[key2].value))
                                else:
                                    print(name2 +' is not in ' + name1)
                            tupkey1.append((name1, dict(tupkey2)))
                        else:
                            tupkey1.append((name1, data1))
                    else:
                        print(name1 + ' is not in ' + name)
                tuples.append((name, dict(tupkey1)))
            else:
                tuples.append((name, data))
        else:
            print(name +' is not in this dataset')

    pprint.pprint('Got GFR dataset: ')
    pprint.pprint(dict(tuples))

    return dict(tuples)

if __name__ == '__main__':

    dataset = read_dcm()


    data = get_gfrdata(dataset)

