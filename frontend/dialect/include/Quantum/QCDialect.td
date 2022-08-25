#ifndef QC_DIALECT
#define QC_DIALECT

include "mlir/IR/OpBase.td"

def QCDialect: Dialect {
    let name = "qcl";
    let summary = "A MLIR dialect for quantum high level languages.";
    let description = [{
        This dialect provides operations specific to quantum high level languages.
    }];
    let cppNamespace = "::mlir::quantum";
}

//===-----------------------------------------------------------------------===//
// Base standalone operation definition.
//===-----------------------------------------------------------------------===//

class QCOp<string mnemonic, list<OpTrait> traits = []> :
        Op<QCDialect, mnemonic, traits>;

#endif // QC_DIALECT



