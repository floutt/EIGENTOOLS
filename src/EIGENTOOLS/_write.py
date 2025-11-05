from math import ceil, nan, isnan
from EIGENTOOLS._read import SNP_Info, Ind_Info
from typing import Literal
import warnings


val_map = {0: 0, 1: 1, 2: 2, nan: 3}


class PackedAncestryMapWriter:
    """
    Object designed for easy writing of a PackedAncestryMap file. This is done on a SNP record basis.

    Attributes:
        _nind (int): Number of individuals based on provided Ind_Info object
        _nsnp (int): Number of variants based on provided SNP_Info object
        _trailingbytes (int): Number of null bytes to write after record
        _isclosed (bool): Boolean indicating if file is closed
        _recordsleft (int): Number of records left to write
    """
    def __init__(self, snp_obj: SNP_Info, ind_obj: Ind_Info,
                 geno_file: (str | None) = None, ind_file: (str | None) = None,
                 snp_file: (str | None) = None,
                 file_prefix: (str | None) = None, write_snp: bool = True,
                 write_ind: bool = True, write_header: bool = True) -> None:
        """
        Initialization method for PackedAncestryMapWriter object.

        Args:
            snp_obj (SNP_Info): Object containing variant information
            ind_obj (Ind_Info): Object containing individual information
            geno_file (str, optional): Output .geno file, cannot be used alongside "file_prefix" parameter
            ind_file (str, optional): Output .ind file, cannot be used alongside "file_prefix" parameter
            snp_file (str, optional): Output .snp file, cannot be used alongside "file_prefix" parameter
            file_prefix (str, optional): Prefix for all PackedAncestryMap files. Will produce a ".ind", ".snp", and ".geno" file. Cannot be used alongside "geno_file", "ind_file", or "snp_file" parameters
            write_snp (boolean): writes variant file if True
            write_ind (boolean): writes individual file if True
            write_header (boolean): writes PackedAncestryMap header if True
        """
        # generate useful metadata
        self._nind = len(ind_obj)
        self._nsnp = len(snp_obj)
        header = ("GENO   %i %i %x %x" % (self._nind, self._nsnp,
                                          ind_obj._hash, snp_obj._hash)).encode()
        min_byte_per_record = ceil(self._nind / 4)
        recordsize = max(min_byte_per_record, len(header))
        self._trailingbytes = bytes(recordsize - min_byte_per_record)
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
            self._fgeno.write(bytes(recordsize - len(header)))

    def write_record(self, dosage_list: list[Literal[0, 1, 2, nan]]) -> None:
        """
        Writes SNP record to PackedAncestryMap

        Args:
            dosage_list (list): List of allelic dosages. Length of list must be equal to the number of individuals.
        """
        shift_by = 6
        record = 0
        # check length
        if len(dosage_list) != self._nind:
            raise ValueError("Length of record should be equal to %i, the number of individuals in the dataset." % self._nind)
        elif self._isclosed:
            raise ValueError("PackedAncestryMapWriter object is closed.")
        for v in dosage_list:
            v = nan if isnan(v) else v  # use proper nan instance
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
        self._fgeno.write(bytes(self._trailingbytes))
        self._recordsleft -= 1

    def close(self) -> None:
        """
        Closes PackedAncestryMap file object and marks object not suitable for writing
        """
        if self._isclosed:
            return None
        else:
            if self._recordsleft != 0:
                warnings.warn("Incomplete number of records written. File may be corrupt.")
            self._fgeno.close()
            self._isclosed = True
