import copy
import os
from math import floor, nan


# hashing functions
def hash_str(s):
    hash_out = 0
    for c in s:
        hash_out *= 23
        hash_out += ord(c)
    return hash_out


# get 32-bit hash from a list of strings.
def hash_list(str_list):
    hash_out = 0
    bit_mask = ((2 ** 32) - 1)
    for s in str_list:
        thash = hash_str(s)
        hash_out *= 17
        hash_out = (hash_out ^ thash) & bit_mask
    return hash_out


# classes used to store the information
class SNP_Info:
    def __init__(self, filename):
        self.var_name = []
        self.chrom = []
        self.pos = []
        self.ref = []
        self.alt = []
        self.cm = []
        self._var_name_to_index = {}
        self._hash = 0

        with open(filename) as f:
            i = 0
            for line in f:
                elems = line.strip().split()
                self.var_name.append(elems[0])
                self.chrom.append(elems[1])
                self.pos.append(int(elems[3]))
                self.ref.append(elems[4])
                self.alt.append(elems[5])
                self.cm.append(float(elems[2]))
                i += 1
        self._var_name_to_index = self._reverse_index()
        self._hash = hash_list(self.var_name)

    def _reverse_index(self):
        out = {}
        for i in range(len(self)):
            out[self.var_name[i]] = i
        return out

    def __getitem__(self, index):
        is_int = isinstance(index, int)
        tmp_obj = copy.copy(self)
        tmp_obj.var_name = [tmp_obj.var_name[index]] if is_int else tmp_obj.var_name[index]
        tmp_obj.chrom = [tmp_obj.chrom[index]] if is_int else tmp_obj.chrom[index]
        tmp_obj.pos = [tmp_obj.pos[index]] if is_int else tmp_obj.pos[index]
        tmp_obj.ref = [tmp_obj.ref[index]] if is_int else tmp_obj.ref[index]
        tmp_obj.alt = [tmp_obj.alt[index]] if is_int else tmp_obj.alt[index]
        tmp_obj.cm = [tmp_obj.cm[index]] if is_int else tmp_obj.cm[index]
        tmp_obj._hash = hash_list(tmp_obj.var_name)
        # redo reverse index
        tmp_obj._var_name_to_index = tmp_obj._reverse_index()
        return tmp_obj

    def __len__(self):
        assert len(self.var_name) == len(self.chrom) == len(self.pos) == len(self.ref) == len(self.alt) == len(self.cm)
        return len(self.var_name)

    def __add__(self, obj2):
        tmp_obj = copy.copy(self)
        tmp_obj.var_name = tmp_obj.var_name + obj2.var_name
        tmp_obj.chrom = tmp_obj.chrom + obj2.chrom
        tmp_obj.pos = tmp_obj.pos + obj2.pos
        tmp_obj.ref = tmp_obj.ref + obj2.ref
        tmp_obj.alt = tmp_obj.alt + obj2.alt
        tmp_obj.cm = tmp_obj.cm + obj2.cm
        tmp_obj._hash = hash_list(tmp_obj.var_name)
        tmp_obj._var_name_to_index = tmp_obj._reverse_index()
        return tmp_obj

    def get_var_name_idx(self, var_name):
        try:
            return self._var_name_to_index[var_name]
        except KeyError:
            raise LookupError("Variant \"" + var_name + "\" not found")

    def write(self, filename):
        with open(filename, "w+") as f:
            for i in range(len(self)):
                f.write("%s\t%s\t%.6f\t%i\t%s\t%s\n" %
                        (self.var_name[i], self.chrom[i], self.cm[i], self.pos[i],
                         self.ref[i], self.alt[i]))


class Ind_Info:
    def __init__(self, filename):
        self.ind_name = []
        self.sex = []
        self.label = []
        self._label_to_idx = {}  # indices for each population
        self._hash = 0

        with open(filename) as f:
            i = 0
            for line in f:
                elems = line.strip().split()
                self.ind_name.append(elems[0])
                self.sex.append(elems[1])
                self.label.append(elems[2])
                i += 1
        self._label_to_idx = self._reverse_index()
        self._hash = hash_list(self.ind_name)

    def __getitem__(self, index):
        is_int = isinstance(index, int)
        tmp_obj = copy.copy(self)
        tmp_obj.ind_name = [tmp_obj.ind_name[index]] if is_int else tmp_obj.ind_name[index]
        tmp_obj.sex = [tmp_obj.sex[index]] if is_int else tmp_obj.sex[index]
        tmp_obj.label = [tmp_obj.label[index]] if is_int else tmp_obj.label[index]
        tmp_obj._hash = hash_list(tmp_obj.ind_name)
        tmp_obj._label_to_idx = tmp_obj._reverse_index()
        return tmp_obj

    def __len__(self):
        assert len(self.ind_name) == len(self.sex) == len(self.label)
        return len(self.ind_name)

    def __add__(self, obj2):
        tmp_obj = copy.copy(self)
        tmp_obj.ind_name = tmp_obj.ind_name + obj2.ind_name
        tmp_obj.sex = tmp_obj.sex + obj2.sex
        tmp_obj.label = tmp_obj.label + obj2.label
        tmp_obj._hash = hash_list(tmp_obj.ind_name)
        tmp_obj._label_to_idx = tmp_obj._reverse_index()
        return tmp_obj

    def _reverse_index(self):
        out = {}
        for i in range(len(self)):
            try:
                out[self.label[i]].append(i)
            except KeyError:
                out[self.label[i]] = [i]
        return out

    def get_label_indices(self, label):
        try:
            return self._label_to_idx[label]
        except KeyError:
            raise LookupError("Label \"" + label + "\" not found.")

    def write(self, filename):
        with open(filename, "w+") as f:
            for i in range(len(self)):
                f.write("%s\t%s\t%s\n" %
                        (self.ind_name[i], self.sex[i], self.label[i]))


class PackedAncestryMap:
    def __init__(self, geno_file=None, ind_file=None, snp_file=None,
                 file_prefix=None, check_hash=True, check_size=True):
        paramter_error_msg = "Inappropriate parametrization. Either only provide a 'file_prefix' parameter or provide parameters for each individual PACKEDANCESTRYMAP file component"
        # parameter handling to allow for init polymorphism
        if file_prefix is None:
            if not ((geno_file is not None) and (ind_file is not None) and
                    (snp_file is not None)):
                raise TypeError(paramter_error_msg)
        else:
            if not ((geno_file is None) and (ind_file is None) and
                    (snp_file is None)):
                raise TypeError(paramter_error_msg)
            else:
                geno_file = file_prefix + ".geno"
                ind_file = file_prefix + ".ind"
                snp_file = file_prefix + ".snp"
        self._GENO_MAP = [0, 1, 2, nan]
        self.snp_info = SNP_Info(snp_file)
        self.ind_info = Ind_Info(ind_file)
        self._fin = open(geno_file, "rb")
        self._filesize = os.path.getsize(geno_file)
        self._recordsize = self._filesize / (len(self.snp_info) + 1)
        assert self._recordsize == floor(self._recordsize)
        self._recordsize = int(self._recordsize)
        self._recordbits = self._recordsize * 8
        raw_header = self._fin.read(self._recordsize).decode()
        # check that header string is in proper format
        if not raw_header.startswith("GENO"):
            raise Exception("Improper .geno filetype")
        self._HEADER = raw_header.replace(chr(0), "")  # remove NULL chars
        # perform optional checks
        header_elems = self._HEADER.split()
        n_ind = int(header_elems[1])
        n_snp = int(header_elems[2])
        ind_hash = int(header_elems[3], 16)
        snp_hash = int(header_elems[4], 16)
        self._i_snp = -1
        self.geno = [0] * len(self.ind_info)
        if (check_hash) & (ind_hash != self.ind_info._hash):
            raise Exception(".ind file hash mismatch. Calculated hash of %x does not match expected hash of %x" % (self.ind_info._hash, ind_hash))
        if (check_hash) & (snp_hash != self.snp_info._hash):
            raise Exception(".snp file hash mismatch. Calculated hash of %x does not match expected hash of %x" % (self.snp_info._hash, snp_hash))
        if (check_size) & (n_ind != len(self.ind_info)):
            raise Exception("Number of individuals in .ind file (n=%i) different from what is expected by .geno file (n=%i)" % (len(self.ind_info), n_ind))
        if (check_size) & (n_snp != len(self.snp_info)):
            raise Exception("Number of SNPs in .snp file (m=%i) different from what is expected by .geno file (m=%i)" % (len(self.snp_info), n_snp))

    def __iter__(self):
        return self

    def get_SNP_Info(self):
        if self._i_snp < 0:
            return None
        return self.snp_info[self._i_snp]

    def _read_record(self):
        rsb = self._recordbits
        snp_record = int.from_bytes(self._fin.read(self._recordsize))
        bit_mask = (2 ** rsb) - 1
        for i in range(len(self.ind_info)):
            dos_tmp = (snp_record << (self._recordbits - rsb)) & bit_mask
            dos_tmp = dos_tmp >> (self._recordbits - 2)
            self.geno[i] = self._GENO_MAP[dos_tmp]
            rsb -= 2

    def __next__(self):
        if self._i_snp == (len(self.snp_info) - 1):
            self._fin.close()
            raise StopIteration
        self._read_record()
        self._i_snp += 1
        return self

    def goto_snp(self, var_name):
        i_var = self.snp_info.get_var_name_idx(var_name)
        self._fin.seek((1 + i_var) * self._recordsize)
        self._read_record()
        self._i_snp = i_var
