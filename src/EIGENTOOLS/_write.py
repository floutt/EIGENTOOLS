from math import ceil, nan
import warnings

val_map = {0: 0, 1: 1, 2: 2, nan: 3}


class PackedAncestryMapWriter:
    def __init__(self, snp_obj, ind_obj, geno_file=None, ind_file=None,
                 snp_file=None, file_prefix=None, write_snp=True,
                 write_ind=True, write_header=True):
        # generate useful metadata
        self._nind = len(ind_obj)
        self._nsnp = len(snp_obj)
        header = ("GENO   %i %i %x %x" % (self._nind, self._nsnp,
                                          ind_obj._hash, snp_obj._hash)).encode()
        self._min_byte_per_record = ceil(self._nind / 4)
        self._recordsize = max(self._min_byte_per_record, len(header))
        self._trailingbytes = bytes(self._recordsize - self._min_byte_per_record)
        self._isclosed = False
        self._recordsleft = self._nsnp

        # write the snp and ind files and the header for the PACKEDANCESTRYMAP
        # parameter handling to allow for init polymorphism
        paramter_error_msg = "Inappropriate parametrization. Either only provide a 'file_prefix' parameter or provide parameters for each individual PACKEDANCESTRYMAP file component"
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

        self._fgeno = open(geno_file, "wb+")

        if write_snp:
            snp_obj.write(snp_file)
        if write_ind:
            ind_obj.write(ind_file)
        # write header
        if write_header:
            self._fgeno.write(header)
            self._fgeno.write(bytes(self._recordsize - len(header)))

    def write_record(self, dosage_list):
        shift_by = 6
        record = 0
        # check length
        if len(dosage_list) != self._nind:
            raise ValueError("Length of record should be equal to %i, the number of individuals in the dataset." % self._nind)
        elif self._isclosed:
            raise ValueError("PackedAncestryMapWriter object is closed.")
        for v in dosage_list:
            try:
                val_map[v]
            except KeyError:
                raise ValueError("Dosages must either be 0,1,2, or nan.")
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
        self._recordsleft -= 1

    def close(self):
        if self._isclosed:
            return None
        else:
            if self._recordsleft != 0:
                warnings.warn("Incomplete number of records written. File may be corrupt.")
            self._fgeno.close()
            self._isclosed = True
