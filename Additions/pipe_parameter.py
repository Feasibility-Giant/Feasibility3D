bl_info = {
    "name": "Pipeline Calculator",
    "author": "Feasibility Giant Company Limited",
    "version": (4, 7, 7),
    "blender": (2, 80, 0),
    "location": "View3D > ToolBar",
    "warning": "",
    "category": "Calc",
}

import bpy
from bpy.props import FloatProperty, BoolProperty
import math
import csv

class Pipeline:
    def __init__(self, length, wall_thickness, fluid_density, fluid_viscosity, atmospheric_pressure):
        # Initialize pipeline properties
        self.length = length
        self.wall_thickness = wall_thickness
        self.fluid_density = fluid_density
        self.fluid_viscosity = fluid_viscosity
        self.atmospheric_pressure = atmospheric_pressure
        self.M_PI = math.pi

    def leak_detection(self, time_interval):
        s_threshold = 0.001
        p_i = self.atmospheric_pressure
        Q = 0

        for t_elapsed in range(0, int(time_interval), 1):
            p = p_i - self.fluid_density * (Q / (self.M_PI * (self.wall_thickness ** 2) / 4) ** 2) / (2 * self.wall_thickness)
            Q_new = math.sqrt((p_i - p) * self.M_PI * (self.wall_thickness ** 2) * 2 / self.fluid_density)

            if abs(Q_new - Q) > s_threshold:
                return True

            p_i = p
            Q = Q_new

        return False

    def leak_localization(self, leak_size, flow_rate):
        d = 2 * leak_size
        v = flow_rate / (self.M_PI * (d ** 2) / 4)
        x = self.fluid_viscosity * leak_size / (self.fluid_density * v * d)
        return x

    def time_of_leak_detection(self, leak_size):
        d = 2 * self.wall_thickness + leak_size
        Q = self.fluid_viscosity * leak_size / (self.fluid_density * d * self.wall_thickness)
        t_detection = self.length / Q
        return t_detection

    def pressure_at_leak_point(self, initial_pressure, leak_size):
        A = self.M_PI * (leak_size ** 2) / 4
        v = math.sqrt(2 * initial_pressure / self.fluid_density)
        p_l = initial_pressure + self.fluid_density * v ** 2 / 2 - self.fluid_density * (
            4 * self.fluid_viscosity * v / (self.fluid_density * A * leak_size)) ** 2 / (
            2 * initial_pressure)
        return p_l

    def leak_volume(self, leak_size):
        d = 2 * self.wall_thickness + leak_size
        Q = self.fluid_viscosity * leak_size / (self.fluid_density * d * self.wall_thickness)
        V = Q * self.length / self.wall_thickness
        return V

class PipelineFlowAssurance:
    def __init__(self, fluid_density, fluid_velocity, pipe_diameter, pipe_length):
        self.fluid_density = fluid_density
        self.fluid_velocity = fluid_velocity
        self.pipe_diameter = pipe_diameter
        self.pipe_length = pipe_length

    def pressure_drop(self, friction_factor):
        return friction_factor * (self.pipe_length / self.pipe_diameter) * (
            self.fluid_density * self.fluid_velocity ** 2) / 2

class CO2Corrosion:
    def __init__(self, corrosion_rate_constant, co2_partial_pressure):
        self.corrosion_rate_constant = corrosion_rate_constant
        self.co2_partial_pressure = co2_partial_pressure

    def corrosion_rate(self):
        return self.corrosion_rate_constant * self.co2_partial_pressure

class H2SCorrosion:
    def __init__(self, corrosion_rate_constant, h2s_partial_pressure):
        self.corrosion_rate_constant = corrosion_rate_constant
        self.h2s_partial_pressure = h2s_partial_pressure

    def corrosion_rate(self):
        return self.corrosion_rate_constant * self.h2s_partial_pressure

class PipelineErosion:
    def __init__(self, fluid_density, fluid_velocity, pipe_diameter, pipe_length):
        self.fluid_density = fluid_density
        self.fluid_velocity = fluid_velocity
        self.pipe_diameter = pipe_diameter
        self.pipe_length = pipe_length

    def erosion_rate(self, erosion_rate_coefficient):
        return erosion_rate_coefficient * (self.fluid_density * self.fluid_velocity ** 3) / (self.pipe_diameter ** 2)

class PipelineHydrates:
    def __init__(self, temperature, pressure):
        self.temperature = temperature
        self.pressure = pressure

    def hydrate_formation_condition(self):
        a = 0.45724
        b = 0.07780
        R = 0.08314
        critical_pressure = 37.5
        critical_temperature = 190.6
        accent_a = (0.42748 * (R ** 2) * (critical_temperature ** 2)) / critical_pressure

        alpha = (1 + (0.37464 + 1.54226 * 0.05134 - 0.26992) * (1 - (self.temperature / critical_temperature)) ** 0.5) ** 2
        a_coef = accent_a * alpha
        b_coef = 0.07780 * R * critical_temperature / critical_pressure

        water_pressure = 0.0314 * math.exp(
            647.3 / self.temperature * (1 - (self.temperature / 647.3)) * (
            -1.2935 + 3.1825 * (self.temperature / 647.3) ** 0.5 - 2.202 * (self.temperature / 647.3) ** 3))
        water_fugacity = water_pressure / self.pressure

        term1 = (self.pressure + b_coef) / (R * self.temperature)
        term2 = (self.pressure + b_coef + 2.414 * a_coef / (R * self.temperature)) / (self.pressure + b_coef - 0.414 * a_coef / (R * self.temperature))

        if term1 <= 0 or term2 <= 0:
            raise ValueError("Logarithm argument must be positive.")

        gas_fugacity = R * self.temperature * math.log(term1) - (a_coef / (2 * math.sqrt(2) * self.temperature)) * math.log(term2)

        if water_fugacity > gas_fugacity:
            return True
        else:
            return False

class PipelineSlugAnalysis:
    def __init__(self, liquid_density, liquid_velocity, pipe_diameter):
        self.liquid_density = liquid_density
        self.liquid_velocity = liquid_velocity
        self.pipe_diameter = pipe_diameter

    def slug_length(self, slug_length_constant):
        return slug_length_constant * math.sqrt(self.liquid_density / self.liquid_velocity)

class PipelineWaxFormation:
    def __init__(self, wax_deposition_rate_constant, wax_appearance_temperature, fluid_temperature):
        self.wax_deposition_rate_constant = wax_deposition_rate_constant
        self.wax_appearance_temperature = wax_appearance_temperature
        self.fluid_temperature = fluid_temperature

    def wax_deposition_rate(self):
        return self.wax_deposition_rate_constant * (self.wax_appearance_temperature - self.fluid_temperature)

class PipeWallShearStress:
    def __init__(self, fluid_density, fluid_velocity, pipe_diameter, pipe_length, roughness):
        self.fluid_density = fluid_density
        self.fluid_velocity = fluid_velocity
        self.pipe_diameter = pipe_diameter
        self.pipe_length = pipe_length
        self.roughness = roughness

    def friction_factor(self, reynolds_number):
        f = 0.005
        while True:
            term1 = (self.roughness / self.pipe_diameter) / (3.7)
            term2 = (2.51 / (reynolds_number * math.sqrt(f)))
            left_side = (1 / math.sqrt(f))
            right_side = (-2 * math.log10(term1 + term2))
            f_new = 1 / (left_side - right_side) ** 2
            if abs(f_new - f) < 0.00001:
                return f_new
            f = f_new

    def wall_shear_stress(self, friction_factor):
        return friction_factor * (self.pipe_length / self.pipe_diameter) * (
            self.fluid_density * self.fluid_velocity ** 2) / 2

class CalculatePipelineParameters(bpy.types.Panel):
    bl_label = "Pipeline Parameters"
    bl_idname = "PT_CalculatePipelineParameters"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'PipeFlow Calculator'
    
    def draw(self, context):
        layout = self.layout
        
        # Fluid Properties
        fluid_section = layout.box()
        fluid_section.label(text="Fluid Properties")
        fluid_section.prop(context.scene, "fluid_density", text="Density (kg/m³)")
        fluid_section.prop(context.scene, "fluid_viscosity", text="Viscosity (Pa.s)")
        
        # Pipe Properties
        pipe_section = layout.box()
        pipe_section.label(text="Pipe Properties")
        pipe_section.prop(context.scene, "pipe_diameter", text="Diameter (m)")
        pipe_section.prop(context.scene, "pipe_length", text="Length (m)")
        pipe_section.prop(context.scene, "pipe_wall_thickness", text="Wall Thickness (m)")
        
        # Leak Detection Properties
        leak_section = layout.box()
        leak_section.label(text="Leak Detection")
        leak_section.prop(context.scene, "leak_size", text="Leak Size (m)")
        leak_section.prop(context.scene, "time_interval", text="Time Interval (s)")
        
        # Leak Localization Properties
        leak_localization_section = layout.box()
        leak_localization_section.label(text="Leak Localization")
        leak_localization_section.prop(context.scene, "flow_rate", text="Flow Rate (m³/s)")
        
        # Erosion Rate Properties
        erosion_section = layout.box()
        erosion_section.label(text="Erosion Rate")
        erosion_section.prop(context.scene, "erosion_rate_coefficient", text="Erosion Rate Coefficient")
        
        # Wax Formation Properties
        wax_section = layout.box()
        wax_section.label(text="Wax Formation")
        wax_section.prop(context.scene, "wax_deposition_rate_constant", text="Deposition Rate Constant")
        wax_section.prop(context.scene, "wax_appearance_temperature", text="Appearance Temperature (°C)")
        wax_section.prop(context.scene, "fluid_temperature", text="Fluid Temperature (°C)")
        
        # Wall Shear Stress Properties
        shear_section = layout.box()
        shear_section.label(text="Wall Shear Stress")
        shear_section.prop(context.scene, "roughness", text="Pipe Roughness (m)")
        
        # Calculate Button
        layout.operator("object.calculate_pipeline_parameters")

class CalculatePipelineParametersOperator(bpy.types.Operator):
    bl_idname = "object.calculate_pipeline_parameters"
    bl_label = "Calculate Pipeline Parameters"
    
    def execute(self, context):
        # Get scene properties
        fluid_density = context.scene.fluid_density
        fluid_viscosity = context.scene.fluid_viscosity
        pipe_diameter = context.scene.pipe_diameter
        pipe_length = context.scene.pipe_length
        pipe_wall_thickness = context.scene.pipe_wall_thickness
        leak_size = context.scene.leak_size
        time_interval = context.scene.time_interval
        flow_rate = context.scene.flow_rate
        erosion_rate_coefficient = context.scene.erosion_rate_coefficient
        wax_deposition_rate_constant = context.scene.wax_deposition_rate_constant
        wax_appearance_temperature = context.scene.wax_appearance_temperature
        fluid_temperature = context.scene.fluid_temperature
        roughness = context.scene.roughness
        
        # Create pipeline objects
        pipeline = Pipeline(pipe_length, pipe_wall_thickness, fluid_density, fluid_viscosity, 101325)
        erosion = PipelineErosion(fluid_density, 0.5, pipe_diameter, pipe_length)  # Assuming fluid_velocity = 0.5 m/s
        wax = PipelineWaxFormation(wax_deposition_rate_constant, wax_appearance_temperature, fluid_temperature)
        shear_stress = PipeWallShearStress(fluid_density, 0.5, pipe_diameter, pipe_length, roughness)  # Assuming fluid_velocity = 0.5 m/s

        # Calculate leak detection
        leak_detected = pipeline.leak_detection(time_interval)
        
        # Calculate erosion rate
        erosion_rate = erosion.erosion_rate(erosion_rate_coefficient)
        
        # Calculate wax deposition rate
        wax_deposition_rate = wax.wax_deposition_rate()
        
        # Calculate wall shear stress
        reynolds_number = (fluid_density * 0.5 * pipe_diameter) / fluid_viscosity  # Assuming fluid_velocity = 0.5 m/s
        friction_factor = shear_stress.friction_factor(reynolds_number)
        wall_shear_stress = shear_stress.wall_shear_stress(friction_factor)
        
        # Display the results in the console
        self.report({'INFO'}, f"Leak Detected: {leak_detected}")
        self.report({'INFO'}, f"Erosion Rate: {erosion_rate:.2f} kg/m²/s")
        self.report({'INFO'}, f"Wax Deposition Rate: {wax_deposition_rate:.2f} kg/m²/s")
        self.report({'INFO'}, f"Wall Shear Stress: {wall_shear_stress:.2f} Pa")
        
        return {'FINISHED'}

def register():
    bpy.utils.register_class(CalculatePipelineParameters)
    bpy.utils.register_class(CalculatePipelineParametersOperator)
    
    bpy.types.Scene.fluid_density = FloatProperty(
        name="Fluid Density",
        description="Density of the fluid in kg/m³",
        default=1000.0
    )
    bpy.types.Scene.fluid_viscosity = FloatProperty(
        name="Fluid Viscosity",
        description="Viscosity of the fluid in Pa.s",
        default=0.001
    )
    bpy.types.Scene.pipe_diameter = FloatProperty(
        name="Pipe Diameter",
        description="Diameter of the pipe in meters",
        default=0.1
    )
    bpy.types.Scene.pipe_length = FloatProperty(
        name="Pipe Length",
        description="Length of the pipe in meters",
        default=100.0
    )
    bpy.types.Scene.pipe_wall_thickness = FloatProperty(
        name="Pipe Wall Thickness",
        description="Wall thickness of the pipe in meters",
        default=0.01
    )
    bpy.types.Scene.leak_size = FloatProperty(
        name="Leak Size",
        description="Size of the leak in meters",
        default=0.01
    )
    bpy.types.Scene.time_interval = FloatProperty(
        name="Time Interval",
        description="Time interval for leak detection in seconds",
        default=10.0
    )
    bpy.types.Scene.flow_rate = FloatProperty(
        name="Flow Rate",
        description="Flow rate for leak localization in cubic meters per second",
        default=0.01
    )
    bpy.types.Scene.erosion_rate_coefficient = FloatProperty(
        name="Erosion Rate Coefficient",
        description="Coefficient for calculating erosion rate",
        default=0.01
    )
    bpy.types.Scene.wax_deposition_rate_constant = FloatProperty(
        name="Wax Deposition Rate Constant",
        description="Constant for calculating wax deposition rate",
        default=0.01
    )
    bpy.types.Scene.wax_appearance_temperature = FloatProperty(
        name="Wax Appearance Temperature",
        description="Temperature at which wax appears in degrees Celsius",
        default=30.0
    )
    bpy.types.Scene.fluid_temperature = FloatProperty(
        name="Fluid Temperature",
        description="Temperature of the fluid in degrees Celsius",
        default=25.0
    )
    bpy.types.Scene.roughness = FloatProperty(
        name="Pipe Roughness",
        description="Roughness of the pipe in meters",
        default=0.0001
    )

def unregister():
    bpy.utils.unregister_class(CalculatePipelineParameters)
    bpy.utils.unregister_class(CalculatePipelineParametersOperator)
    
    del bpy.types.Scene.fluid_density
    del bpy.types.Scene.fluid_viscosity
    del bpy.types.Scene.pipe_diameter
    del bpy.types.Scene.pipe_length
    del bpy.types.Scene.pipe_wall_thickness
    del bpy.types.Scene.leak_size
    del bpy.types.Scene.time_interval
    del bpy.types.Scene.flow_rate
    del bpy.types.Scene.erosion_rate_coefficient
    del bpy.types.Scene.wax_deposition_rate_constant
    del bpy.types.Scene.wax_appearance_temperature
    del bpy.types.Scene.fluid_temperature
    del bpy.types.Scene.roughness

if __name__ == "__main__":
    register()