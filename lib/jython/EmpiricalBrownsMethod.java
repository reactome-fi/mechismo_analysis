import org.python.core.*;

public class EmpiricalBrownsMethod extends java.lang.Object {
    static String[] jpy$mainProperties = new String[] {"python.modules.builtin", "exceptions:org.python.core.exceptions"};
    static String[] jpy$proxyProperties = new String[] {"python.modules.builtin", "exceptions:org.python.core.exceptions", "python.options.showJavaExceptions", "true"};
    static String[] jpy$packages = new String[] {};
    
    public static class _PyInner extends PyFunctionTable implements PyRunnable {
        private static PyObject s$0;
        private static PyObject i$1;
        private static PyObject i$2;
        private static PyObject f$3;
        private static PyObject i$4;
        private static PyObject f$5;
        private static PyObject f$6;
        private static PyObject f$7;
        private static PyObject f$8;
        private static PyObject f$9;
        private static PyObject i$10;
        private static PyObject s$11;
        private static PyFunctionTable funcTable;
        private static PyCode c$0_EmpiricalBrownsMethod;
        private static PyCode c$1_lambda;
        private static PyCode c$2_TransformData;
        private static PyCode c$3_CalculateCovariances;
        private static PyCode c$4_CombinePValues;
        private static PyCode c$5_KostsMethod;
        private static PyCode c$6_KostPolyFit;
        private static PyCode c$7_CalculateKostCovariance;
        private static PyCode c$8_main;
        private static void initConstants() {
            s$0 = Py.newString("\012Created on Sun Jun 21 12:22:12 2015\012\012@author:  William Poole: wpoole@caltech.edu\012");
            i$1 = Py.newInteger(2);
            i$2 = Py.newInteger(0);
            f$3 = Py.newFloat(2.0);
            i$4 = Py.newInteger(1);
            f$5 = Py.newFloat(4.0);
            f$6 = Py.newFloat(1.0);
            f$7 = Py.newFloat(3.263);
            f$8 = Py.newFloat(0.71);
            f$9 = Py.newFloat(0.027);
            i$10 = Py.newInteger(3);
            s$11 = Py.newString("/Users/joshuaburkhart/SoftwareProjects/CombiningDependentPvaluesUsingEBM/Python/EmpiricalBrownsMethod.py");
            funcTable = new _PyInner();
            c$0_EmpiricalBrownsMethod = Py.newCode(3, new String[] {"data_matrix", "p_values", "extra_info", "covar_matrix"}, "/Users/joshuaburkhart/SoftwareProjects/CombiningDependentPvaluesUsingEBM/Python/EmpiricalBrownsMethod.py", "EmpiricalBrownsMethod", false, false, funcTable, 4, null, null, 0, 17);
            c$1_lambda = Py.newCode(1, new String[] {"x"}, "/Users/joshuaburkhart/SoftwareProjects/CombiningDependentPvaluesUsingEBM/Python/EmpiricalBrownsMethod.py", "<lambda>", false, false, funcTable, 5, null, new String[] {"s"}, 0, 17);
            c$2_TransformData = Py.newCode(1, new String[] {"data_vector", "W", "m", "_[2]", "d", "sd", "x", "_[1]", "s"}, "/Users/joshuaburkhart/SoftwareProjects/CombiningDependentPvaluesUsingEBM/Python/EmpiricalBrownsMethod.py", "TransformData", false, false, funcTable, 6, new String[] {"s"}, null, 1, 17);
            c$3_CalculateCovariances = Py.newCode(1, new String[] {"data_matrix", "covar_matrix", "_[1]", "transformed_data_matrix", "f"}, "/Users/joshuaburkhart/SoftwareProjects/CombiningDependentPvaluesUsingEBM/Python/EmpiricalBrownsMethod.py", "CalculateCovariances", false, false, funcTable, 7, null, null, 0, 17);
            c$4_CombinePValues = Py.newCode(3, new String[] {"covar_matrix", "p_values", "extra_info", "Var", "_[1]", "p_brown", "Expected", "x", "df_brown", "p", "p_fisher", "m", "j", "i", "df_fisher", "cov_sum", "c"}, "/Users/joshuaburkhart/SoftwareProjects/CombiningDependentPvaluesUsingEBM/Python/EmpiricalBrownsMethod.py", "CombinePValues", false, false, funcTable, 8, null, null, 0, 17);
            c$5_KostsMethod = Py.newCode(3, new String[] {"data_matrix", "p_values", "extra_info", "covar_matrix"}, "/Users/joshuaburkhart/SoftwareProjects/CombiningDependentPvaluesUsingEBM/Python/EmpiricalBrownsMethod.py", "KostsMethod", false, false, funcTable, 9, null, null, 0, 17);
            c$6_KostPolyFit = Py.newCode(1, new String[] {"cor", "a2", "a1", "a3"}, "/Users/joshuaburkhart/SoftwareProjects/CombiningDependentPvaluesUsingEBM/Python/EmpiricalBrownsMethod.py", "KostPolyFit", false, false, funcTable, 10, null, null, 0, 17);
            c$7_CalculateKostCovariance = Py.newCode(1, new String[] {"data_matrix", "m", "covar", "j", "i", "p_val", "covar_matrix", "cor"}, "/Users/joshuaburkhart/SoftwareProjects/CombiningDependentPvaluesUsingEBM/Python/EmpiricalBrownsMethod.py", "CalculateKostCovariance", false, false, funcTable, 11, null, null, 0, 17);
            c$8_main = Py.newCode(0, new String[] {}, "/Users/joshuaburkhart/SoftwareProjects/CombiningDependentPvaluesUsingEBM/Python/EmpiricalBrownsMethod.py", "main", false, false, funcTable, 12, null, null, 0, 16);
        }
        
        
        public PyCode getMain() {
            if (c$8_main == null) _PyInner.initConstants();
            return c$8_main;
        }
        
        public PyObject call_function(int index, PyFrame frame) {
            switch (index){
                case 0:
                return _PyInner.__listcomprehension$1(frame);
                case 1:
                return _PyInner.__listcomprehension$2(frame);
                case 2:
                return _PyInner.__listcomprehension$3(frame);
                case 3:
                return _PyInner.__listcomprehension$4(frame);
                case 4:
                return _PyInner.EmpiricalBrownsMethod$5(frame);
                case 5:
                return _PyInner.lambda$6(frame);
                case 6:
                return _PyInner.TransformData$7(frame);
                case 7:
                return _PyInner.CalculateCovariances$8(frame);
                case 8:
                return _PyInner.CombinePValues$9(frame);
                case 9:
                return _PyInner.KostsMethod$10(frame);
                case 10:
                return _PyInner.KostPolyFit$11(frame);
                case 11:
                return _PyInner.CalculateKostCovariance$12(frame);
                case 12:
                return _PyInner.main$13(frame);
                default:
                return null;
            }
        }
        
        private static PyObject __listcomprehension$1(PyFrame frame) {
            // Temporary Variables
            PyObject t$0$PyObject, t$1$PyObject, t$2$PyObject, t$3$PyObject;
            
            t$0$PyObject = new PyList(new PyObject[] {});
            frame.setlocal(7, t$0$PyObject.__getattr__("append"));
            t$2$PyObject = frame.getlocal(0).__iter__();
            while ((t$3$PyObject = t$2$PyObject.__iternext__()) != null) {
                frame.setlocal(4, t$3$PyObject);
                frame.getlocal(7).__call__(frame.getlocal(4)._sub(frame.getlocal(2))._div(frame.getlocal(5)));
            }
            return t$0$PyObject;
        }
        
        private static PyObject __listcomprehension$2(PyFrame frame) {
            // Temporary Variables
            PyObject t$0$PyObject, t$1$PyObject, t$2$PyObject, t$3$PyObject;
            
            t$0$PyObject = new PyList(new PyObject[] {});
            frame.setlocal(3, t$0$PyObject.__getattr__("append"));
            t$2$PyObject = frame.getderef(0).__iter__();
            while ((t$3$PyObject = t$2$PyObject.__iternext__()) != null) {
                frame.setlocal(6, t$3$PyObject);
                frame.getlocal(3).__call__(frame.getlocal(1).__call__(frame.getlocal(6)));
            }
            return t$0$PyObject;
        }
        
        private static PyObject __listcomprehension$3(PyFrame frame) {
            // Temporary Variables
            PyObject t$0$PyObject, t$1$PyObject, t$2$PyObject, t$3$PyObject;
            
            t$0$PyObject = new PyList(new PyObject[] {});
            frame.setlocal(2, t$0$PyObject.__getattr__("append"));
            t$2$PyObject = frame.getlocal(0).__iter__();
            while ((t$3$PyObject = t$2$PyObject.__iternext__()) != null) {
                frame.setlocal(4, t$3$PyObject);
                frame.getlocal(2).__call__(frame.getglobal("TransformData").__call__(frame.getlocal(4)));
            }
            return t$0$PyObject;
        }
        
        private static PyObject __listcomprehension$4(PyFrame frame) {
            // Temporary Variables
            PyObject t$0$PyObject, t$1$PyObject, t$2$PyObject, t$3$PyObject;
            
            t$0$PyObject = new PyList(new PyObject[] {});
            frame.setlocal(4, t$0$PyObject.__getattr__("append"));
            t$2$PyObject = frame.getlocal(1).__iter__();
            while ((t$3$PyObject = t$2$PyObject.__iternext__()) != null) {
                frame.setlocal(9, t$3$PyObject);
                frame.getlocal(4).__call__(frame.getglobal("np").__getattr__("log").__call__(frame.getlocal(9)).__neg__());
            }
            return t$0$PyObject;
        }
        
        private static PyObject EmpiricalBrownsMethod$5(PyFrame frame) {
            frame.setlocal(3, frame.getglobal("CalculateCovariances").__call__(frame.getlocal(0)));
            return frame.getglobal("CombinePValues").__call__(frame.getlocal(3), frame.getlocal(1), frame.getlocal(2));
        }
        
        private static PyObject lambda$6(PyFrame frame) {
            return i$1.__neg__()._mul(frame.getglobal("np").__getattr__("log").__call__(frame.getglobal("ECDF").__call__(frame.getderef(0)).__call__(frame.getlocal(0))));
        }
        
        private static PyObject TransformData$7(PyFrame frame) {
            frame.setlocal(2, frame.getglobal("np").__getattr__("mean").__call__(frame.getlocal(0)));
            frame.setlocal(5, frame.getglobal("np").__getattr__("std").__call__(frame.getlocal(0)));
            frame.setderef(0, __listcomprehension$1(frame));
            frame.setlocal(1, new PyFunction(frame.f_globals, new PyObject[] {}, c$1_lambda, new PyObject[] {frame.getclosure(0)}));
            return frame.getglobal("np").__getattr__("array").__call__(__listcomprehension$2(frame));
        }
        
        private static PyObject CalculateCovariances$8(PyFrame frame) {
            frame.setlocal(3, frame.getglobal("np").__getattr__("array").__call__(__listcomprehension$3(frame)));
            frame.setlocal(1, frame.getglobal("np").__getattr__("cov").__call__(frame.getlocal(3)));
            return frame.getlocal(1);
        }
        
        private static PyObject CombinePValues$9(PyFrame frame) {
            // Temporary Variables
            PyObject t$0$PyObject, t$1$PyObject, t$2$PyObject, t$3$PyObject;
            
            // Code
            frame.setlocal(11, frame.getglobal("int").__call__(frame.getlocal(0).__getattr__("shape").__getitem__(i$2)));
            frame.setlocal(14, f$3._mul(frame.getlocal(11)));
            frame.setlocal(6, f$3._mul(frame.getlocal(11)));
            frame.setlocal(15, i$2);
            t$0$PyObject = frame.getglobal("range").__call__(frame.getlocal(11)).__iter__();
            while ((t$1$PyObject = t$0$PyObject.__iternext__()) != null) {
                frame.setlocal(13, t$1$PyObject);
                t$2$PyObject = frame.getglobal("range").__call__(frame.getlocal(13)._add(i$4), frame.getlocal(11)).__iter__();
                while ((t$3$PyObject = t$2$PyObject.__iternext__()) != null) {
                    frame.setlocal(12, t$3$PyObject);
                    frame.setlocal(15, frame.getlocal(15).__iadd__(frame.getlocal(0).__getitem__(new PyTuple(new PyObject[] {frame.getlocal(13), frame.getlocal(12)}))));
                }
            }
            frame.setlocal(3, f$5._mul(frame.getlocal(11))._add(i$1._mul(frame.getlocal(15))));
            frame.setlocal(16, frame.getlocal(3)._div(f$3._mul(frame.getlocal(6))));
            frame.setlocal(8, f$3._mul(frame.getlocal(6)._pow(i$1))._div(frame.getlocal(3)));
            if (frame.getlocal(8)._gt(frame.getlocal(14)).__nonzero__()) {
                frame.setlocal(8, frame.getlocal(14));
                frame.setlocal(16, f$6);
            }
            frame.setlocal(7, f$3._mul(frame.getglobal("sum").__call__(__listcomprehension$4(frame))));
            frame.setlocal(5, frame.getglobal("chi2_cdf").__call__(frame.getlocal(8), f$6._mul(frame.getlocal(7))._div(frame.getlocal(16))));
            frame.setlocal(10, frame.getglobal("chi2_cdf").__call__(frame.getlocal(14), f$6._mul(frame.getlocal(7))));
            if (frame.getlocal(2).__nonzero__()) {
                return new PyTuple(new PyObject[] {frame.getlocal(5), frame.getlocal(10), frame.getlocal(16), frame.getlocal(8)});
            }
            else {
                return frame.getlocal(5);
            }
        }
        
        private static PyObject KostsMethod$10(PyFrame frame) {
            frame.setlocal(3, frame.getglobal("CalculateKostCovariance").__call__(frame.getlocal(0)));
            return frame.getglobal("CombinePValues").__call__(new PyObject[] {frame.getlocal(3), frame.getlocal(1), frame.getlocal(2)}, new String[] {"extra_info"});
        }
        
        private static PyObject KostPolyFit$11(PyFrame frame) {
            // Temporary Variables
            PyObject[] t$0$PyObject__;
            
            // Code
            t$0$PyObject__ = org.python.core.Py.unpackSequence(new PyTuple(new PyObject[] {f$7, f$8, f$9}), 3);
            frame.setlocal(2, t$0$PyObject__[0]);
            frame.setlocal(1, t$0$PyObject__[1]);
            frame.setlocal(3, t$0$PyObject__[2]);
            return frame.getlocal(2)._mul(frame.getlocal(0))._add(frame.getlocal(1)._mul(frame.getlocal(0)._pow(i$1)))._add(frame.getlocal(3)._mul(frame.getlocal(0)._pow(i$10)));
        }
        
        private static PyObject CalculateKostCovariance$12(PyFrame frame) {
            // Temporary Variables
            PyObject[] t$0$PyObject__;
            PyObject t$0$PyObject, t$1$PyObject, t$2$PyObject, t$3$PyObject;
            
            // Code
            frame.setlocal(1, frame.getlocal(0).__getattr__("shape").__getitem__(i$2));
            frame.setlocal(6, frame.getglobal("np").__getattr__("zeros").__call__(new PyTuple(new PyObject[] {frame.getlocal(1), frame.getlocal(1)})));
            t$0$PyObject = frame.getglobal("range").__call__(frame.getlocal(1)).__iter__();
            while ((t$1$PyObject = t$0$PyObject.__iternext__()) != null) {
                frame.setlocal(4, t$1$PyObject);
                t$2$PyObject = frame.getglobal("range").__call__(frame.getlocal(4)._add(i$4), frame.getlocal(1)).__iter__();
                while ((t$3$PyObject = t$2$PyObject.__iternext__()) != null) {
                    frame.setlocal(3, t$3$PyObject);
                    t$0$PyObject__ = org.python.core.Py.unpackSequence(frame.getglobal("pearsonr").__call__(frame.getlocal(0).__getitem__(new PyTuple(new PyObject[] {frame.getlocal(4), new PySlice(null, null, null)})), frame.getlocal(0).__getitem__(new PyTuple(new PyObject[] {frame.getlocal(3), new PySlice(null, null, null)}))), 2);
                    frame.setlocal(7, t$0$PyObject__[0]);
                    frame.setlocal(5, t$0$PyObject__[1]);
                    frame.setlocal(2, frame.getglobal("KostPolyFit").__call__(frame.getlocal(7)));
                    frame.getlocal(6).__setitem__(new PyTuple(new PyObject[] {frame.getlocal(4), frame.getlocal(3)}), frame.getlocal(2));
                    frame.getlocal(6).__setitem__(new PyTuple(new PyObject[] {frame.getlocal(3), frame.getlocal(4)}), frame.getlocal(2));
                }
            }
            return frame.getlocal(6);
        }
        
        private static PyObject main$13(PyFrame frame) {
            frame.setglobal("__file__", s$11);
            
            // Temporary Variables
            PyObject[] t$0$PyObject__;
            
            // Code
            /* 
            Created on Sun Jun 21 12:22:12 2015
            
            @author:  William Poole: wpoole@caltech.edu
             */
            frame.setlocal("np", org.python.core.imp.importOneAs("numpy", frame));
            t$0$PyObject__ = org.python.core.imp.importFrom("statsmodels.distributions.empirical_distribution", new String[] {"ECDF"}, frame);
            frame.setlocal("ECDF", t$0$PyObject__[0]);
            t$0$PyObject__ = null;
            t$0$PyObject__ = org.python.core.imp.importFrom("scipy.special", new String[] {"chdtrc"}, frame);
            frame.setlocal("chi2_cdf", t$0$PyObject__[0]);
            t$0$PyObject__ = null;
            t$0$PyObject__ = org.python.core.imp.importFrom("scipy.stats", new String[] {"pearsonr"}, frame);
            frame.setlocal("pearsonr", t$0$PyObject__[0]);
            t$0$PyObject__ = null;
            frame.setlocal("EmpiricalBrownsMethod", new PyFunction(frame.f_globals, new PyObject[] {frame.getname("False")}, c$0_EmpiricalBrownsMethod));
            frame.setlocal("TransformData", new PyFunction(frame.f_globals, new PyObject[] {}, c$2_TransformData));
            frame.setlocal("CalculateCovariances", new PyFunction(frame.f_globals, new PyObject[] {}, c$3_CalculateCovariances));
            frame.setlocal("CombinePValues", new PyFunction(frame.f_globals, new PyObject[] {frame.getname("False")}, c$4_CombinePValues));
            frame.setlocal("KostsMethod", new PyFunction(frame.f_globals, new PyObject[] {frame.getname("False")}, c$5_KostsMethod));
            frame.setlocal("KostPolyFit", new PyFunction(frame.f_globals, new PyObject[] {}, c$6_KostPolyFit));
            frame.setlocal("CalculateKostCovariance", new PyFunction(frame.f_globals, new PyObject[] {}, c$7_CalculateKostCovariance));
            return Py.None;
        }
        
    }
    public static void moduleDictInit(PyObject dict) {
        dict.__setitem__("__name__", new PyString("EmpiricalBrownsMethod"));
        Py.runCode(new _PyInner().getMain(), dict, dict);
    }
    
    public static void main(String[] args) throws java.lang.Exception {
        String[] newargs = new String[args.length+1];
        newargs[0] = "EmpiricalBrownsMethod";
        java.lang.System.arraycopy(args, 0, newargs, 1, args.length);
        Py.runMain(EmpiricalBrownsMethod._PyInner.class, newargs, EmpiricalBrownsMethod.jpy$packages, EmpiricalBrownsMethod.jpy$mainProperties, null, new String[] {"EmpiricalBrownsMethod"});
    }
    
}
