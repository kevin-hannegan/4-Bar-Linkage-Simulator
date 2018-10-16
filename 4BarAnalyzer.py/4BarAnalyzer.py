from math import sin, cos, asin, acos, sqrt, pow, pi, radians, degrees
from time import sleep
import tkinter as tk
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

class InputArea(tk.Frame):
    def __init__(self, master):
        self.master = master
        tk.Frame.__init__(self, master)
        self.frame = tk.Frame(self.master, bg='white')
        vcmd = (self.master.register(self.validateNumeric), '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        self.l1Label = tk.Label(self.frame, text="Rocker/Crank Link Length (L1)", justify="left", bg='white')
        self.l1Label.grid(row=0, column=0, sticky="w")
        self.l1Value = tk.StringVar()
        self.l1Entry = tk.Entry(self.frame, validate="key", validatecommand=vcmd, textvariable=self.l1Value)
        self.l1Entry.grid(row=0, column=1, sticky="w")
        self.l2Label = tk.Label(self.frame, text="Coupler Link Length (L2)", justify="left", bg='white')
        self.l2Label.grid(row=1, column=0, sticky="w")
        self.l2Value = tk.StringVar()
        self.l2Entry = tk.Entry(self.frame, validate="key", validatecommand=vcmd, textvariable=self.l2Value)
        self.l2Entry.grid(row=1, column=1, sticky="w")
        self.l3Label = tk.Label(self.frame, text="Follower Link Length (L3)", justify="left", bg='white')
        self.l3Label.grid(row=2, column=0, sticky="w")
        self.l3Value = tk.StringVar()
        self.l3Entry = tk.Entry(self.frame, validate="key", validatecommand=vcmd, textvariable=self.l3Value)
        self.l3Entry.grid(row=2, column=1, sticky="w")
        self.l4Label = tk.Label(self.frame, text="Ground/Stationary Link Length (L4)", justify="left", bg='white')
        self.l4Label.grid(row=3, column=0, sticky="w")
        self.l4Value = tk.StringVar()
        self.l4Entry = tk.Entry(self.frame, validate="key", validatecommand=vcmd, textvariable=self.l4Value)
        self.l4Entry.grid(row=3, column=1, sticky="w")
        self.theta4Label = tk.Label(self.frame, text="Ground/Stationary Link Angle (Theta_4, degrees)", justify="left", bg='white')
        self.theta4Label.grid(row=4, column=0, sticky="w")
        self.theta4Value = tk.StringVar()
        self.theta4Value.set(0)
        self.theta4Entry = tk.Entry(self.frame, validate="key", validatecommand=vcmd, textvariable=self.theta4Value)
        self.theta4Entry.grid(row=4, column=1, sticky="w")
        self.trackerPointLabel = tk.Label(self.frame, text="Tracker Point Coordinates (Origin bisects coupler)", justify="left", bg='white')
        self.trackerPointLabel.grid(row=5, column=0, sticky="w")
        self.trackerPointXValue = tk.StringVar()
        self.trackerPointYValue = tk.StringVar()
        self.trackerPointXValue.set(0)
        self.trackerPointYValue.set(0)
        self.trackerPointFrame = tk.Frame(self.frame)
        self.trackerPointXEntry = tk.Entry(self.trackerPointFrame, validate="key", validatecommand=vcmd, textvariable=self.trackerPointXValue, width=7)
        self.trackerPointYEntry = tk.Entry(self.trackerPointFrame, validate="key", validatecommand=vcmd, textvariable=self.trackerPointYValue, width=7)
        self.trackerPointXLabel = tk.Label(self.trackerPointFrame, text="X:", bg='white')
        self.trackerPointYLabel = tk.Label(self.trackerPointFrame, text="Y:", bg='white')
        self.trackerPointXLabel.grid(row=0, column=0)
        self.trackerPointXEntry.grid(row=0, column=1)
        self.trackerPointYLabel.grid(row=0, column=2)
        self.trackerPointYEntry.grid(row=0, column=3)
        self.trackerPointFrame.grid(row=5, column=1, sticky="w")
        self.upperBoundLabel = tk.Label(self.frame, text="Theta_1 Upper Bound Angle (Degrees)", justify="left", bg='white')
        self.upperBoundValue = tk.StringVar()
        self.upperBoundValue.set("0")
        self.upperBoundEntry = tk.Entry(self.frame, validate="key", validatecommand=vcmd, textvariable=self.upperBoundValue)
        self.upperBoundLabel.grid(row=6, column=0, sticky="w")
        self.upperBoundEntry.grid(row=6, column=1, sticky="w")
        self.lowerBoundLabel = tk.Label(self.frame, text="Theta_1 Lower Bound Angle (Degrees)", justify="left", bg='white')
        self.lowerBoundValue = tk.StringVar()
        self.lowerBoundValue.set("0")
        self.lowerBoundEntry = tk.Entry(self.frame, validate="key", validatecommand=vcmd, textvariable=self.lowerBoundValue)
        self.lowerBoundLabel.grid(row=7, column=0, sticky="w")
        self.lowerBoundEntry.grid(row=7, column=1, sticky="w")
        self.startAngleLabel = tk.Label(self.frame, text="Simulation Start Angle (Degrees)", justify="left", bg='white')
        self.startAngleValue = tk.StringVar()
        self.startAngleValue.set("90")
        self.startAngleEntry = tk.Entry(self.frame, validate="key", validatecommand=vcmd, textvariable=self.startAngleValue)
        self.startAngleLabel.grid(row=8, column=0, sticky="w")
        self.startAngleEntry.grid(row=8, column=1, sticky="w")
        self.incrementLabel = tk.Label(self.frame, text="Simulation Step Increment (Degrees)", justify="left", bg='white')
        self.incrementValue = tk.StringVar()
        self.incrementValue.set("5")
        self.incrementEntry = tk.Entry(self.frame, validate="key", validatecommand=vcmd, textvariable=self.incrementValue)
        self.incrementLabel.grid(row=9, column=0, sticky="w")
        self.incrementEntry.grid(row=9, column=1, sticky="w")
        self.messageText = tk.StringVar()
        self.messageBox = tk.LabelFrame(self.frame, text="Messages", padx=5, pady=20, bg='white')
        self.messageText.set("Program Initializing...")
        self.messageLabel = tk.Label(self.messageBox, textvariable=self.messageText, anchor="center", justify="left", bg='white')
        self.messageLabel.grid(row=0, column=0, sticky="news")
        self.messageBox.grid(row=10, column=0, columnspan=2, sticky="news")
        self.dataBox = tk.LabelFrame(self.frame, text="Data", padx=5, pady=20, bg='white')
        self.theta1Degrees = tk.StringVar()
        self.theta1LabelText = tk.Label(self.dataBox, text="Theta_1:", justify="left", bg='white')
        self.theta1CalcText = tk.Label(self.dataBox, textvariable=self.theta1Degrees, bg='white')
        self.theta2Degrees = tk.StringVar()
        self.theta2LabelText = tk.Label(self.dataBox, text="Theta_2:", justify="left", bg='white')
        self.theta2CalcText = tk.Label(self.dataBox, textvariable=self.theta2Degrees, bg='white')
        self.theta3Degrees = tk.StringVar()
        self.theta3LabelText = tk.Label(self.dataBox, text="Theta_3:", justify="left", bg='white')
        self.theta3CalcText = tk.Label(self.dataBox, textvariable=self.theta3Degrees, bg='white')
        self.theta1LabelText.grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.theta1CalcText.grid(row=0, column=1, padx=5, pady=5, sticky="w")
        self.theta2LabelText.grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.theta2CalcText.grid(row=1, column=1, padx=5, pady=5, sticky="w")
        self.theta3LabelText.grid(row=2, column=0, padx=5, pady=5, sticky="w")
        self.theta3CalcText.grid(row=2, column=1, padx=5, pady=5, sticky="w")
        self.dataBox.grid(row=11, column=0, columnspan=2, sticky="news")
        self.frame.grid(row=0, column=1, sticky="N", padx=5, pady=5)

    def validateNumeric(self, action, index, value_if_allowed, prior_value, text, validation_type, trigger_type, widget_name):
        if text == '':
            return True
        if text in '0123456789.-+' or text == '':
            try:
                float(value_if_allowed)
                return True
            except ValueError:
                return False
        else:
            return False



class PlotArea(tk.Frame):
    def __init__(self, master):
        self.master = master
        tk.Frame.__init__(self, master)
        self.figure = Figure(figsize=(5,5), dpi=100)
        self.axes = self.figure.add_subplot(111)
        self.axes.spines['right'].set_visible(False)
        self.axes.spines['top'].set_visible(False)
        self.canvas = FigureCanvasTkAgg(self.figure, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0)
        self.canvas._tkcanvas.grid(row=0, column=0)

class ButtonArea(tk.Frame):
    def __init__(self, master):
        self.master = master
        tk.Frame.__init__(self, master)
        self.frame = tk.Frame(self.master)
        self.simulationToggleText = tk.StringVar()
        self.simulationToggleText.set("Press to Run Simulation")
        self.simulationToggleButton = tk.Button(self.frame, textvariable=self.simulationToggleText, bg="green", fg="white", command=self.toggleSimulation, width=20)
        self.simulationOnSwitch = tk.BooleanVar() # use trace on this variable from parent to track state and impact simulation
        self.simulationOnSwitch.set(False)
        self.exitButton = tk.Button(self.frame, text="Exit", command=self.killProgram, width=20)
        self.killBool = tk.BooleanVar()
        self.killBool.set(False)
        self.frame.grid(row=1, column=1)
        self.simulationToggleButton.grid(row=0, column=0)
        self.exitButton.grid(row=0, column=1)

    def toggleSimulation(self):
        if self.simulationOnSwitch.get() == False:
            self.simulationToggleButton.configure(bg="red")
            self.simulationToggleText.set("Pause Simulation")
            self.simulationOnSwitch.set(True)
        else:
            self.simulationToggleButton.configure(bg="green")
            self.simulationToggleText.set("Press to Run Simulation")
            self.simulationOnSwitch.set(False)

    def killProgram(self):
        self.killBool.set(True)

class AnalysisArea(tk.Frame):
    def __init__(self, master):
        self.master = master
        tk.Frame.__init__(self, master)
        self.figure = Figure(figsize=(10, 3), dpi=100)
        self.theta2Axes = self.figure.add_subplot(211)
        self.theta3Axes = self.figure.add_subplot(212, sharex=self.theta2Axes)
        self.theta2Line, = self.theta2Axes.plot([], [], 'r-')
        self.theta3Line, = self.theta3Axes.plot([], [], 'b-')
        self.figure.suptitle('theta_2, theta_3 as function of crank angle')
        self.theta2Axes.set_ylabel('theta_2\n (degrees)', rotation=0, labelpad=30)
        self.theta3Axes.set_ylabel(ylabel='theta_3\n (degrees)', rotation=0, labelpad=30)
        self.theta3Axes.set_xlabel(xlabel='theta_1 (degrees)')
        self.figure.subplots_adjust(hspace=0)
        self.canvas = FigureCanvasTkAgg(self.figure, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=2, column=0)
        self.canvas._tkcanvas.grid(row=2, column=0, columnspan=2)

class Log:
    def __init__(self):
        self.theta1 = []
        self.theta2 = []
        self.theta3 = []
        self.trackerPath = []
        self.n = [] # number of completed iterations

    def logData(self, theta1, theta2, theta3, trackerPoint, n):
        self.theta1.append(theta1)
        self.theta2.append(theta2)
        self.theta3.append(theta3)
        self.trackerPath.append(trackerPoint.get())
        self.n.append(n)

    def printLog(self):
        print("n\t theta1\t theta2\t theta3\t t1-t2\t")
        for i in range(0, len(self.n)):
            print("{}\t {}\t {}\t {}\t {}\t".format(self.n[i], self.theta1[i], self.theta2[i], self.theta3[i], self.theta1[i] - self.theta2[i]))



class Linkage:
    def __init__(self, L1=.8, L2=1, L3=1, L4=1, CouplerOffsetX=0, CouplerOffsetY=0):
        self.L1 = L1 # Length of crank linkage
        self.L2 = L2 # Length of coupler linkage
        self.L3 = L3 # length of follower linkage
        self.L4 = L4 # length of ground linkage
        self.couplerOffsetX = CouplerOffsetX
        self.couplerOffsetY = CouplerOffsetY
        self.theta1 = pi / 2 # use radians, right angle to start location
        self.precision = .000001
        self.theta4 = 0 # gets set by entry
        self.theta3 = 0 # initial guess
        self.theta2 = 0 # initial guess
        self.origin = Point(0,0)
        self.anchorPoint = Point(self.L4 * cos(self.theta4), self.L4 * sin(self.theta4))
        self.n = 0
        self.log = Log()
        self.validJoints = False
        self.trackerPoint = Point()
        self.J1 = Point()
        self.J2 = Point()
        

        
        

    def get_theta1(self):
        return self.theta1

    def set_theta1(self, theta1):
        theta1Old = self.theta1 # get theta1 from previous step, probably obsolete
        self.theta1 = theta1 # sets theta1
        self.solveGeometry()

    def solveGeometry(self): 
        self.theta2 = self.newton(self.theta2)
        self.theta3 = self.findTheta3(self.theta2)
        self.findJoints()
        self.log.logData(degrees(self.theta1), degrees(self.theta2), degrees(self.theta3), self.trackerPoint, self.n)
        self.n += 1

    def newton(self, theta2Guess, iterationLevel=0):
        convergenceError = self.precision + 10 # setting convergenceError > precision so loop executes, dummy value
        previousGuess = theta2Guess
        # for debugging
        if iterationLevel > 0:
            print("Iteration Level: {}".format(iterationLevel))
            sleep(.01)
        try:
            counter = 0
            while abs(convergenceError) > self.precision:
                f = self.fourBarPhysics(theta2Guess)
                df = self.dFourBarPhysics(theta2Guess, self.theta1)
                theta2Guess = theta2Guess - f / df
                convergenceError = f
                #print("newton counter: {}, Convergence Error: {}".format(counter, convergenceError))
                counter += 1
                if counter > 50:
                    raise ValueError("Newton's method not converging")
        except ValueError:
            # after 10 (arbitrary number) tries where newton's method fails to converge due to discontinuous functions
            # Try an  runge-kutta method approximation
            if iterationLevel < 10 :
                
                self.newton(previousGuess - 2 * pi / 3600, iterationLevel + 1) # adds a tenth degree and tries again, test retry system
            else:
                if self.n > 2:
                    theta2Guess = self.rungeKutta(self.dFourBarPhysics, previousGuess, self.theta1)
        return theta2Guess

    def rungeKutta(self, funcDTheta, previousValue, theta1):
        h  = self.log.theta1[len(self.log.theta1)-1] - self.log.theta1[len(self.log.theta1)-2] # step size across theata1, copy from previous step size
        k1 = h * funcDTheta(previousValue, theta1)
        k2 = h * funcDTheta(previousValue + h/2, theta1 + k1/2)
        k3 = h * funcDTheta(previousValue + h/2, theta1 + k2/2)
        k4 = h * funcDTheta(previousValue + h, theta1 + k3)
        thetaGuess = previousValue + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)
        return thetaGuess


    def fourBarPhysics(self, theta2):
        argument = self.L1 * cos(self.theta1) / self.L3 + self.L2 * cos(self.theta1 - theta2) / self.L3 - self.L4 * cos(self.theta4) / self.L3
        #print("Acos domain test - theta2: {}, argument {}".format(theta2, argument))
        intermediateCalc = acos(self.L1 * cos(self.theta1) / self.L3 + self.L2 * cos(self.theta1 - theta2) / self.L3 - self.L4 * cos(self.theta4) / self.L3)
        root = self.L1 * sin(self.theta1) + self.L2 * sin(self.theta1 - theta2) - self.L3 * sin(intermediateCalc) - self.L4 * sin(self.theta4)
        return root
    
    def dFourBarPhysics(self, theta2, theta1):
        try:
            numerator = (self.L2 / self.L3) * (sin(theta1 - theta2)) * (self.L1 * cos(theta1) + self.L2 * cos(theta1 - theta2) - self.L4 * cos(self.theta4))
            denominator = sqrt(1 - pow((self.L1 * cos(theta1) / self.L3 + self.L2 * cos(theta1 - theta2) / self.L3 - self.L4 * cos(self.theta4) / self.L3) ,2))
            offset = - self.L2 * cos(theta1 - theta2)
            dTheta2 = numerator / denominator + offset
        except ValueError:
            # near discontinuities, approximate with backward difference method, error correction to occur on future steps when newton's method finds root
            # TODO: lots of coding, but may be able to central difference method to decrease error, but requires stepping acrros discontinuity
            stepsize = radians(self.log.theta1[len(self.log.theta1) - 2]) - radians(self.log.theta1[len(self.log.theta1) - 1])
            rise = radians(self.log.theta1[len(self.log.theta2) - 2]) - radians(self.log.theta1[len(self.log.theta2) - 1])
            dTheta2 = rise / stepsize
        return dTheta2

    def findTheta3(self, theta2):
        try:
            # analytical solution to solive to theta 3
            argument = (self.L1 * cos(self.theta1) + self.L2 * cos(self.theta1 - theta2) - self.L4 * cos(self.theta4)) / self.L3
            #print("theta3 argument: {}".format(argument))
            theta3 = acos(argument)
        except ValueError:
            # use runge-kutta t 
            theta3 = self.rungeKutta(self.dTheta3, theta2, self.theta1)
        return theta3

    def dTheta3(self, theta2, theta1):
        try:
            # analytical derivatie of theta3
            numerator = self.L1 * sin(theta1) / self.L3 + self.L2 * sin(theta1 - theta2) / self.L3
            denominator = sqrt(1 - pow((1 / self.L3) * (self.L1 * cos(theta1) + self.L2 * cos(theta1 - theta2) - self.L4 * cos(self.theta4) ), 2))
            dTheta3 = numerator / denominator
        except ValueError:
            # use backwards difference method to approximate 
            stepsize = radians(self.log.theta1[len(self.log.theta1) - 2]) - radians(self.log.theta1[len(self.log.theta1) - 1])
            rise = radians(self.log.theta1[len(self.log.theta3) - 2]) - radians(self.log.theta1[len(self.log.theta3) - 1])
            dTheta3 = rise / stepsize
        return dTheta3


    def findJoints(self):
        self.anchorPoint = Point(self.L4 * cos(self.theta4), self.L4 * sin(self.theta4))
        self.J1 = Point(self.L1 * cos(self.theta1), self.L1 * sin(self.theta1))
        self.J2 = Point(self.L3 * cos(self.theta3) + self.L4 * cos(self.theta4), self.L3 * sin(self.theta3) + self.L4 * sin(self.theta4))
        xPrime = self.L1 * cos(self.theta1) + self.L2 * cos(self.theta1 - self.theta2)
        yPrime = self.L1 * sin(self.theta1) + self.L2 * sin(self.theta1 - self.theta2)
        self.J2Prime = Point(xPrime, yPrime)
        self.isJointValid(self.J2, self.J2Prime)
        self.trackerOrigin = self.bisectPoints(self.J1, self.J2)
        trackerPointX = self.trackerOrigin.x - self.couplerOffsetY * sin(self.theta1 - self.theta2) + self.couplerOffsetX * cos(self.theta1 - self.theta2)
        trackerPointY = self.trackerOrigin.y + self.couplerOffsetY * cos(self.theta1 - self.theta2) + self.couplerOffsetX * sin(self.theta1 - self.theta2)
        self.trackerPoint = Point(trackerPointX, trackerPointY)

    def bisectPoints(self, p1, p2):
        xMid = (p1.x + p2.x) / 2
        yMid = (p1.y + p2.y) / 2
        return Point(xMid, yMid)


    def isJointValid(self, p1, p2):
        distance = sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2))
        if distance < sqrt(self.precision):
            self.validJoints = True
            return True
        else:
            print("Joint displacement error: {}".format(distance))
            self.validJoints = False
            return False

    def draw(self):
        self.rocker = self.drawLine(self.origin, self.J1)
        self.coupler = self.drawLine(self.J1, self.J2)
        self.follower = self.drawLine(self.J2, self.anchorPoint)
        self.offset = self.drawLine(self.trackerOrigin, self.trackerPoint)
        self.offset.set_marker(None)
        return self.rocker, self.coupler, self.follower, self.offset

    def drawLine(self, p1, p2):
        xdata, ydata = self.getLineData(p1, p2)
        line = plt.Line2D(xdata, ydata, lw=4, marker='.', markersize=10, color='black', zorder=0)
        return line
    
    def drawTracker(self):
        xdata, ydata = self.getTrackerData()
        trackerLine = plt.Line2D(xdata, ydata, linestyle='dashed', color='red', linewidth=2)
        return trackerLine
    
    def getTrackerData(self):
    
        xdata = [i[0] for i in self.log.trackerPath]
        ydata = [i[1] for i in self.log.trackerPath]
        return xdata, ydata
       


    def getLineData(self, p1, p2):
        xdata = (p1.x, p2.x)
        ydata = (p1.y, p2.y)
        return xdata, ydata
    
    def isGrashof(self):
        nonRockers = [self.L2, self.L3, self.L4]
        longestLink = max(nonRockers)
        remainingLinkSum = sum(nonRockers) - longestLink
        rockerPlusLongest = longestLink + self.L1
        if rockerPlusLongest < remainingLinkSum:
            return True
        return False

class Point:
    def __init__(self, x=0 , y=0):
        self.x = x
        self.y = y

    def get(self):
        return [self.x, self.y]

    def set(self, x , y):
        self.x = x
        self.y = y



class Pin:
    def __init__(self, anchorPoint, size, orientation):
        self.anchorPoint = anchorPoint
        self.size = size
        self.orientation = radians(orientation) # convert from degrees
    
    def drawPin(self):
        pinHole = plt.Circle((self.anchorPoint.x, self.anchorPoint.y), radius=self.size/5, fc='black')
        mount = plt.Circle((self.anchorPoint.x, self.anchorPoint.y), radius=self.size, ec='black')
        baseBracket = [[self.size * -1, 0], [-2* self.size,-2 *self.size], [self.size * 2, -2 * self.size], [self.size, 0]]
        bracketPoints = []
        for x,y in baseBracket:
            #rotate around origin
            x_prime = x * cos(self.orientation) - y * sin(self.orientation)
            y_prime = y * cos(self.orientation) + x * sin(self.orientation)
            # transform to anchor location
            x_prime = x_prime + self.anchorPoint.x
            y_prime = y_prime + self.anchorPoint.y
            newPoint = [x_prime, y_prime]
            bracketPoints.append(newPoint)
        bracket = plt.Polygon(bracketPoints, closed=None, fill=None, edgecolor="black")
        basePt1 = [bracketPoints[1][0] + self.size / 2 * sin(self.orientation), bracketPoints[1][1] - self.size / 2 * cos(self.orientation)]
        basePt2 = [bracketPoints[2][0] + self.size / 2 * sin(self.orientation), bracketPoints[2][1] - self.size / 2 * cos(self.orientation)]
        basePoints = [bracketPoints[1], basePt1, basePt2, bracketPoints[2]]
        base = plt.Polygon(basePoints, closed=True, fill="black")
        return pinHole, mount, bracket, base







class Simulator(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        # crate Gui geometry
        self.configure(bg='white')
        self.container = tk.Frame(self, bg='white', padx=20, pady=20)
        self.container.grid(row=0, column=0)
        self.plotArea = PlotArea(self.container)
        self.inputArea = InputArea(self.container)
        self.buttonArea = ButtonArea(self.container)

        self.plotArea.grid(row=0, column=0, rowspan=2)
        self.inputArea.grid(row=0, column=0)
        self.buttonArea.grid(row=1, column=1)
        # self.analysisArea = AnalysisArea(self.container) add back for debugging
        # self.analysisArea.grid(row=2, column=0, columnspan=2)

        # Initialize Linkage
        self.linkage = Linkage()
        self.updateInputs()
        #self.updateLinkage() # not needed in initialization, call in run loop/animation
        self.setBounds()
        self.addOrigin()
        self.addAnchor()  
        self.linkage.solveGeometry()
        self.tracker = self.linkage.drawTracker()
        self.plotArea.axes.add_line(self.tracker)

        # Add event handlers
        self.inputArea.l1Value.trace("w", self.changeGeometryCallback)
        self.inputArea.l2Value.trace("w", self.changeGeometryCallback)
        self.inputArea.l3Value.trace("w", self.changeGeometryCallback)
        self.inputArea.l4Value.trace("w", self.changeGeometryCallback)
        self.inputArea.theta4Value.trace("w", self.changeGeometryCallback)
        self.buttonArea.simulationOnSwitch.trace("w", self.run)
        self.buttonArea.killBool.trace("w", self.exit)
        self.inputArea.upperBoundValue.trace("w", self.changeParameterCallback)
        self.inputArea.lowerBoundValue.trace("w", self.changeParameterCallback)
        self.inputArea.incrementValue.trace("w", self.changeParameterCallback)
        self.inputArea.startAngleValue.trace("w", self.changeParameterCallback)
        self.inputArea.trackerPointXValue.trace("w", self.changeGeometryCallback)
        self.inputArea.trackerPointYValue.trace("w", self.changeGeometryCallback)
        
 
        # draw linkage and update GUI

        self.drawLinkage()
        self.updateMessage()
        self.changeParameterCallback()
        self.changeGeometryCallback()
        self.direction = 1
        # initialize analysis

        #self.plotAnalysis()

        # add tracking arc
        
        

    def updateMessage(self):
        if self.linkage.isGrashof():
            self.inputArea.messageText.set("\nGrashof condition satisfied.\nRocker will act as 360 degree crank.\nSet each boundary angle to 0 for continuous 360 Degree Motion.\n")
        else:
            self.inputArea.messageText.set("\nGrashof condition failed.\n\n Simulation will fail when input angle is not probably bounded.\n")

    def incrementTheta1(self):
        if abs(self.upperBound) > 0.000001 and abs(self.lowerBound) > 0.000001:
            if abs(self.linkage.theta1 - radians(self.upperBound)) <= radians(self.increment) or abs(self.linkage.theta1 - radians(self.lowerBound)) <= radians(self.increment):
                self.direction *= -1
        increment = self.increment * self.direction
        self.linkage.set_theta1(self.linkage.theta1 + radians(increment))

    def updateLinkage(self):
        self.linkage.L1 = float(self.inputArea.l1Value.get())
        self.linkage.L2 = float(self.inputArea.l2Value.get())
        self.linkage.L3 = float(self.inputArea.l3Value.get())
        self.linkage.L4 = float(self.inputArea.l4Value.get())
        self.linkage.theta4 = radians(float(self.inputArea.theta4Value.get()))
        self.linkage.couplerOffsetX = float(self.inputArea.trackerPointXValue.get())
        self.linkage.couplerOffsetY = float(self.inputArea.trackerPointYValue.get())


    def updateInputs(self):
        self.inputArea.l1Value.set(self.linkage.L1)
        self.inputArea.l2Value.set(self.linkage.L2)
        self.inputArea.l3Value.set(self.linkage.L3)
        self.inputArea.l4Value.set(self.linkage.L4)
        self.inputArea.theta4Value.set(self.linkage.theta4)

    def setBounds(self):
        extents = [self.linkage.origin, self.linkage.anchorPoint, self.linkage.J1, self.linkage.J2, self.linkage.trackerPoint]
        ymin, ymax = self.plotArea.axes.get_ylim()
        xmin, xmax = self.plotArea.axes.get_xlim()
        for point in extents:
            if point.x * 1.1 > xmax:
                xmax = max(xmax + (xmax - xmin) / 10, point.x * 1.1)
            if point.x * 1.1 < xmin:
                xmin = min(xmin - (xmax + xmin) / 10, point.x * 1.1)
            if point.y * 1.1 > ymax:
                ymax = max(ymax + (ymax - ymin) / 10, point.y * 1.1)
            if point.y < ymin:
                ymin = min(ymin - (ymax + ymin) / 10, point.y * 1.1)
        #print("xmin: {}, xmax: {}, ymin: {}, ymax: {}".format(xmin, xmax, ymin, ymax)) # for debugging only
        self.plotArea.axes.set_ylim([ymin, ymax])
        self.plotArea.axes.set_xlim([xmin, xmax])
        self.plotArea.axes.set_aspect('equal')

    def addOrigin(self):
        self.originPin = Pin(self.linkage.origin, self.linkage.L4 / 20, 0)
        oPinHole, oMount, oBracket, oBase = self.originPin.drawPin()
        self.originMount = self.plotArea.axes.add_patch(oMount)
        self.originBracket = self.plotArea.axes.add_patch(oBracket)
        self.originBase = self.plotArea.axes.add_patch(oBase)
        self.originPinHole = self.plotArea.axes.add_patch(oPinHole)
        return self.originPin

    def addAnchor(self):
        self.anchorPin = Pin(self.linkage.anchorPoint, self.linkage.L4 / 20, 0)
        oPinHole, oMount, oBracket, oBase = self.anchorPin.drawPin()
        self.anchorMount = self.plotArea.axes.add_patch(oMount)
        self.anchorBracket = self.plotArea.axes.add_patch(oBracket)
        self.anchorBase = self.plotArea.axes.add_patch(oBase)
        self.anchorPinHole = self.plotArea.axes.add_patch(oPinHole)
        return self.anchorPin

    def redrawAnchor(self):
        self.anchorMount.remove()
        self.anchorBracket.remove()
        self.anchorBase.remove()
        self.anchorPinHole.remove()
        self.anchorPin = self.addAnchor()
        self.originBase.remove()
        self.originBracket.remove()
        self.originPinHole.remove()
        self.originMount.remove()
        self.originPin = self.addOrigin()

    def drawLinkage(self):
        self.rocker, self.coupler, self.follower, self.offset = self.linkage.draw()
        self.plotArea.axes.add_line(self.rocker)
        self.plotArea.axes.add_line(self.coupler)
        self.plotArea.axes.add_line(self.follower)
        self.plotArea.axes.add_line(self.offset)
        self.trackerPoint = self.plotArea.axes.scatter(self.linkage.trackerPoint.x, self.linkage.trackerPoint.y, marker='.', color='red', s=200, zorder=10)




    def redrawLinkage(self):

        self.redrawAnchor()
        rockerData = self.linkage.getLineData(self.linkage.origin, self.linkage.J1)
        self.rocker.set_data(rockerData)
        couplerData = self.linkage.getLineData(self.linkage.J1, self.linkage.J2)
        self.coupler.set_data(couplerData)
        followerData = self.linkage.getLineData(self.linkage.J2, self.linkage.anchorPoint)
        self.follower.set_data(followerData)
        offsetData = self.linkage.getLineData(self.linkage.trackerOrigin, self.linkage.trackerPoint)
        self.offset.set_data(offsetData)
        self.trackerPoint.set_offsets([self.linkage.trackerPoint.x, self.linkage.trackerPoint.y])
        if len(self.linkage.log.trackerPath) * self.increment <= 360 + self.increment:
            self.redrawTracker()
            self.setBounds()
        self.plotArea.canvas.draw()
        self.plotArea.canvas.flush_events()

    def redrawTracker(self):
        trackerData = self.linkage.getTrackerData()
        self.tracker.set_data(trackerData)

    def clearTracker(self):
        self.linkage.log.trackerPath = [self.linkage.trackerPoint.get()]

    def changeGeometryCallback(self, *args):
        self.updateLinkage()
        if self.linkage.isGrashof():
            self.linkage.solveGeometry()
            self.redrawLinkage()
            if len(self.linkage.log.trackerPath) > 0:
                self.clearTracker()
        self.updateMessage()

    def changeParameterCallback(self, *args):
        self.upperBound = float(self.inputArea.upperBoundValue.get())
        self.lowerBound = float(self.inputArea.lowerBoundValue.get())
        self.startAngle = float(self.inputArea.startAngleValue.get())
        self.increment = float(self.inputArea.incrementValue.get())


    def updateDataDisplay(self):
        self.inputArea.theta1Degrees.set(str(round(degrees(self.linkage.theta1)  % 360 ,1)))
        self.inputArea.theta2Degrees.set(str(round(degrees(self.linkage.theta2)  % 360, 1)))
        self.inputArea.theta3Degrees.set(str(round(degrees(self.linkage.theta3)  % 360, 1)))

    def run(self, *args):
        self.updateLinkage()
        while self.buttonArea.simulationOnSwitch.get() == True:
            self.incrementTheta1()

            if self.linkage.validJoints == False:
                if self.linkage.isGrashof():
                    self.inputArea.messageText.set("\n\nSimulation Failure.\n\nReduce increment angle size and Try again.\n\nLikely failure of Both newton's method and Runge-Kutta to produce accurate estimtates")
                else:
                    self.inputArea.messageText.set("\n\nSimulation Failure.\n\nCheck Grashof condition and crank angle bounds.")
                self.buttonArea.simulationOnSwitch.set(False)
                return

            self.redrawLinkage()
            # self.plotAnalysis() # for debugging 
            self.updateDataDisplay()


    """ Only using for debugging geometry """
    def plotAnalysis(self):
        self.analysisArea.theta2Line.set_ydata(self.linkage.log.theta2)
        self.analysisArea.theta2Line.set_xdata(self.linkage.log.theta1)
        self.analysisArea.theta3Line.set_ydata(self.linkage.log.theta2)
        self.analysisArea.theta3Line.set_xdata(self.linkage.log.theta1)
        self.analysisArea.canvas.draw()
        self.analysisArea.canvas.flush_events()
            

    def exit(self, *args):
        self.destroy()


if __name__ == '__main__':
    app = Simulator()
    app.mainloop()

