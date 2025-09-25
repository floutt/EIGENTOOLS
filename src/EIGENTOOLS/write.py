from math import ceil, nan

val_map = {0: 0, 1: 1, 2: 2, nan: 3}


class PackedAncestryMapWriter:
    def __init__(self, file_prefix, snp_obj, ind_obj):
        # generate useful metadata
        self._nind = len(snp_obj)
        self._nsnp = len(ind_obj)
        header = ("GENO   %i %i %x %x" % (self._nind, self._nsnp,
                                          ind_obj._hash, snp_obj._hash)).encode()
        self._min_byte_per_record = ceil(self._nind / 4)
        self._recordsize = max(self._min_byte_per_record, len(header))
        self._trailingbytes = bytes(self._recordsize - self._min_byte_per_record)

        # write the snp and ind files and the header for the PACKEDANCESTRYMAP
        snp_file = file_prefix + ".snp"
        ind_file = file_prefix + ".ind"
        geno_file = file_prefix + ".geno"
        self._fgeno = open(geno_file, "wb+")

        snp_obj.write(snp_file)
        ind_obj.write(ind_file)
        # write header
        self._fgeno.write(header)
        self._fgeno.write(bytes(self._recordsize - len(header)))

    def write_record(self, dosage_list):
        shift_by = 6
        record = 0
        for v in dosage_list:
            if (shift_by == 0):
                record = record + (val_map[v] << shift_by)
                self._fgeno.write(bytes([record]))
                record = 0
                shift_by = 6
                continue
            record = record + (val_map[v] << shift_by)
            shift_by -= 2
        if shift_by != 6:
            self._fgeno.write(bytes([record]))
        # write trailing bytes
        self._fgeno.write(bytes(self._recordsize - self._min_byte_per_record))
