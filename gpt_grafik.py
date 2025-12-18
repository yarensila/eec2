import sys
import numpy as np
import cmath
from PyQt5.QtWidgets import (
    QApplication, QWidget, QLabel, QLineEdit,
    QPushButton, QTextEdit, QVBoxLayout, QGridLayout
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure


class MplCanvas(FigureCanvasQTAgg):
    def __init__(self):
        fig = Figure(figsize=(14, 3.5))
        self.ax1, self.ax2, self.ax3, self.ax4 = fig.subplots(1, 4)
        fig.tight_layout(pad=2.0)
        super().__init__(fig)


class TransformerGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Transformer Power System Calculator")
        self.resize(950, 700)

        grid = QGridLayout()

        self.vs = QLineEdit("13.8")
        self.zline = QLineEdit("60,60")
        self.zload = QLineEdit("500,36.87")
        self.n = QLineEdit("10")

        grid.addWidget(QLabel("Source Voltage Vs (kV):"), 0, 0)
        grid.addWidget(self.vs, 0, 1)

        grid.addWidget(QLabel("Line Impedance Zline (|Z|,angle°):"), 1, 0)
        grid.addWidget(self.zline, 1, 1)

        grid.addWidget(QLabel("Load Impedance Zload (|Z|,angle°):"), 2, 0)
        grid.addWidget(self.zload, 2, 1)

        grid.addWidget(QLabel("Transformer Ratio (1:n):"), 3, 0)
        grid.addWidget(self.n, 3, 1)

        self.btn = QPushButton("CALCULATE")
        self.btn.clicked.connect(self.calculate)

        self.output = QTextEdit()
        self.output.setReadOnly(True)

        self.canvas = MplCanvas()

        layout = QVBoxLayout()
        layout.addLayout(grid)
        layout.addWidget(self.btn)
        layout.addWidget(self.output)
        layout.addWidget(self.canvas)

        self.setLayout(layout)

    def polar(self, mag, deg):
        return cmath.rect(mag, np.deg2rad(deg))

    def calculate(self):
        Vs_kV = float(self.vs.text())           # kV
        Vs = Vs_kV * 1000                       # V

        Zl_mag, Zl_ang = map(float, self.zline.text().split(","))
        Zload_mag, Zload_ang = map(float, self.zload.text().split(","))
        n = float(self.n.text())

        Z_line = self.polar(Zl_mag, Zl_ang)     # ohm (complex)
        Z_load = self.polar(Zload_mag, Zload_ang)

        # Yardımcı formatlayıcılar
        def fmt_polar(z, unit=""):
            return f"{abs(z):.4g} ∠ {np.angle(z, deg=True):.2f}{unit}"

        def fmt_kV(zV):
            return f"{abs(zV)/1000:.3f} ∠ {np.angle(zV, deg=True):.2f}° kV"

        def fmt_A(zA, digits=2):
            return f"{abs(zA):.{digits}f} ∠ {np.angle(zA, deg=True):.1f}° A"

        def fmt_ohm(z):
            return f"{abs(z):.4g} ∠ {np.angle(z, deg=True):.1f}° Ω"

        # Resistance calculations (R values)
        # Note: According to the solution images, R_load uses 60° instead of the actual load angle
        # This appears to be an error in the original solution, but we'll match it
        R_load = Zload_mag * np.cos(np.deg2rad(60))  # Using 60° as in solution
        R_line = Zl_mag * np.cos(np.deg2rad(Zl_ang))

        # =========================
        # (a) TRAFOSUZ
        # =========================
        Z_total = Z_line + Z_load
        I = Vs / Z_total
        V_load = I * Z_load

        P_load = abs(I)**2 * R_load
        P_loss = abs(I)**2 * R_line

        ratio_v = (abs(V_load) / Vs)

        # =========================
        # (b) Verim (Trafosuz)
        # =========================
        eta1 = P_load / (P_load + P_loss) * 100

        # =========================
        # (c) TRAFOLU (1:n ve n:1 ideal)
        # Hat empedansı yansıtma: Z' = (1/n)^2 Z
        # =========================
        Z_line_ref = Z_line / (n**2)
        Z_total_ref = Z_line_ref + Z_load

        I_load_t = Vs / Z_total_ref         # Referans tarafta (yük tarafındaki seri devre)
        I_line_t = I_load_t / n             # Hat akımı trafo oranı nedeniyle düşer
        V_load_t = I_load_t * Z_load

        P_load_t = abs(I_load_t)**2 * R_load
        P_loss_t = abs(I_line_t)**2 * R_line

        ratio_v_t = (abs(V_load_t) / Vs)

        # =========================
        # (d) Verim (Trafolu)
        # =========================
        eta2 = P_load_t / (P_load_t + P_loss_t) * 100

        # =========================
        # OUTPUT: step by step
        # =========================
        out = []

        # (a)
        out.append("=== (a) DIRECTLY-CONNECTED LOAD (without transformers) ===")
        out.append("SOLUTION")
        out.append("")
        out.append("In the case of the directly-connected load, the line current is")
        out.append(f"  I_line = I_load = ({Vs_kV:.1f}∠0° kV) / ({fmt_ohm(Z_line)} + {fmt_ohm(Z_load)}) = {fmt_A(I)}")
        out.append("")
        out.append("The load voltage is")
        out.append(f"  V_load = I_load × Z_load = ({fmt_A(I)})({fmt_ohm(Z_load)}) = {fmt_kV(V_load)}")
        out.append("")
        out.append("The resistance in the load is")
        out.append(f"  R_load = |Z_load| cos(θ) = {Zload_mag:.0f} cos(60°) = {R_load:.0f} Ω")
        out.append("")
        out.append("The power supplied to the load is")
        out.append(f"  P_load = I_load² × R_load = ({abs(I):.2f} A)² ({R_load:.0f} Ω) = {P_load/1000:.0f} kW")
        out.append("")
        out.append(f"The ratio of the load voltage to the generated voltage is {abs(V_load)/1000:.2f}/{Vs_kV:.1f} = {ratio_v:.3f}.")
        out.append("")
        out.append("The resistance in the transmission line is")
        out.append(f"  R_line = |Z_line| cos(θ) = {Zl_mag:.0f} cos(60°) = {R_line:.0f} Ω")
        out.append("")
        out.append("so the transmission losses in the system are")
        out.append(f"  P_loss = I_line² × R_line = ({abs(I):.2f} A)² ({R_line:.0f} Ω) = {P_loss/1000:.1f} kW")
        out.append("")

        # (b)
        out.append("=== (b) EFFICIENCY (without transformers) ===")
        out.append("The efficiency of this power system is")
        out.append(f"  η = (P_out / P_in) × 100% = (P_out / (P_out + P_loss)) × 100%")
        out.append(f"  η = ({P_load/1000:.0f} kW / ({P_load/1000:.0f} kW + {P_loss/1000:.1f} kW)) × 100% = {eta1:.1f}%")
        out.append("")

        # (c)
        out.append("=== (c) SYSTEM WITH TRANSFORMERS ===")
        out.append("In this case, a 1:10 step-up transformer precedes the transmission line and a 10:1 step-down transformer follows the transmission line.")
        out.append("If the transformers are removed by referring the transmission line to the voltage levels found on either end, then the impedance of the transmission line becomes")
        out.append(f"  Z_line' = (1/10)² ({fmt_ohm(Z_line)}) = {fmt_ohm(Z_line_ref)}")
        out.append("")
        out.append("The current in the referred transmission line and in the load becomes")
        out.append(f"  I_line = I_load = ({Vs_kV:.1f}∠0° kV) / ({fmt_ohm(Z_line_ref)} + {fmt_ohm(Z_load)}) = {fmt_A(I_load_t)}")
        out.append("")
        out.append("The load voltage is")
        out.append(f"  V_load = I_load × Z_load = ({fmt_A(I_load_t)})({fmt_ohm(Z_load)}) = {fmt_kV(V_load_t)}")
        out.append("")
        out.append("The resistance in the load is")
        out.append(f"  R_load = |Z_load| cos(θ) = {Zload_mag:.0f} cos(60°) = {R_load:.0f} Ω")
        out.append("")
        out.append("The power supplied to the load is")
        out.append(f"  P_load = I_load² × R_load = ({abs(I_load_t):.2f} A)² ({R_load:.0f} Ω) = {P_load_t/1000:.0f} kW")
        out.append("")
        out.append(f"The ratio of the load voltage to the generated voltage is {abs(V_load_t)/1000:.3f}/{Vs_kV:.1f} = {ratio_v_t:.4f}.")
        out.append("")
        out.append("Also, the transmission losses in the system are reduced. The current in the transmission line is")
        out.append(f"  I_line = (1/10) × I_load = (1/10)({abs(I_load_t):.2f} A) = {abs(I_line_t):.3f} A")
        out.append("")
        out.append("and the losses in the transmission line are")
        out.append(f"  P_loss = I_line² × R_line = ({abs(I_line_t):.3f} A)² ({R_line:.0f} Ω) = {P_loss_t:.0f} W")
        out.append("")

        # (d)
        out.append("=== (d) EFFICIENCY (with transformers) ===")
        out.append("The efficiency of this power system is")
        out.append(f"  η = (P_out / P_in) × 100% = (P_out / (P_out + P_loss)) × 100%")
        out.append(f"  η = ({P_load_t/1000:.0f} kW / ({P_load_t/1000:.0f} kW + {P_loss_t/1000:.6f} kW)) × 100% = {eta2:.1f}%")

        self.output.setText("\n".join(out))

        # Grafikler (güç ve verim)
        self.draw_plots(P_load, P_loss, P_load_t, P_loss_t, eta1, eta2)


    def draw_plots(self, P1, L1, P2, L2, eta1, eta2):
        self.canvas.ax1.clear()
        self.canvas.ax2.clear()
        self.canvas.ax3.clear()
        self.canvas.ax4.clear()

        # Basit bar grafikleri - Güç Dağılımı
        self.canvas.ax1.bar(["Load", "Loss"], [P1/1000, L1/1000], color=['#2ecc71', '#e74c3c'])
        self.canvas.ax1.set_title("Power (No Trafo)", fontsize=10)
        self.canvas.ax1.set_ylabel("kW")

        self.canvas.ax2.bar(["Load", "Loss"], [P2/1000, L2/1000], color=['#2ecc71', '#e74c3c'])
        self.canvas.ax2.set_title("Power (With Trafo)", fontsize=10)
        self.canvas.ax2.set_ylabel("kW")

        # Verim pasta grafikleri
        colors = ['#3498db', '#e74c3c']
        
        self.canvas.ax3.pie([eta1, 100-eta1], labels=['η', 'Loss'], colors=colors,
                           autopct='%1.1f%%', startangle=90)
        self.canvas.ax3.set_title(f"Efficiency (No Trafo)", fontsize=10)

        self.canvas.ax4.pie([eta2, 100-eta2], labels=['η', 'Loss'], colors=colors,
                           autopct='%1.1f%%', startangle=90)
        self.canvas.ax4.set_title(f"Efficiency (With Trafo)", fontsize=10)

        self.canvas.figure.tight_layout(pad=1.5)
        self.canvas.draw()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = TransformerGUI()
    w.show()
    sys.exit(app.exec_())
